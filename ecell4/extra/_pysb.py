import numbers, operator, itertools, functools
from collections.abc import Iterable

import ecell4_base.core
import pysb

class PySBModel(object):

    def __init__(self, name, _export=False):
        self._export = _export
        self.__model = pysb.Model(name=name, _export=self._export)

    def add_component(self, cls, *args, **kwargs):
        obj = cls(*args, _export=self._export, **kwargs)
        self.__model.add_component(obj)
        return obj

    def add_monomer(self, *args, **kwargs):
        return self.add_component(pysb.Monomer, *args, **kwargs)

    def add_parameter(self, *args, **kwargs):
        return self.add_component(pysb.Parameter, *args, **kwargs)

    def __add_rule(self, *args, **kwargs):
        return self.add_component(pysb.Rule, *args, **kwargs)

    def add_rule(self, name, left, right, rate_forward, rate_reverse=None, *args, **kwargs):
        if len(left) == 0:
            left = None
        else:
            left = sum(left[1: ], left[0])
        if len(right) == 0:
            right = None
        else:
            right = sum(right[1: ], right[0])
        if rate_reverse is None:
            rule = operator.rshift(left, right)
        else:
            rule = operator.or_(left, right)
        return self.__add_rule(name, rule, rate_forward, rate_reverse, *args, **kwargs)

    def add_observable(self, *args, **kwargs):
        return self.add_component(pysb.Observable, *args, **kwargs)

    def add_initial(self, *args, **kwargs):
        obj = pysb.Initial(*args, _export=self._export, **kwargs)
        self.__model.add_initial(obj)
        return obj

    @property
    def monomers(self):
        return self.__model.monomers

    @property
    def parameters(self):
        return self.__model.parameters

    @property
    def observables(self):
        return self.__model.observables

    def get(self):
        #XXX copy?
        return self.__model

    def run(self, t):
        if isinstance(t, numbers.Real):
            import numpy
            t = numpy.linspace(0, t, 100)
        import pysb.simulator
        return pysb.simulator.ScipyOdeSimulator(self.__model, t).run()

    def plot(self, *, t=None, result=None):
        assert t is not None or result is not None
        result = result or self.run(t)

        import matplotlib.pylab as plt
        fig = plt.figure()
        for obj in result._model.observables:
            plt.plot(result.tout[0], result.observables[obj.name], label=obj.name)
        plt.xlabel('Time')
        plt.ylabel('Amount')
        plt.legend(loc='best')
        plt.show()

def list_monomers(obj, monomers=None):
    monomers = monomers or {}
    if isinstance(obj, (ecell4_base.core.NetworkModel, ecell4_base.core.NetfreeModel)):
        for rr in obj.reaction_rules():
            monomers = list_monomers(rr, monomers)
        return monomers
    elif isinstance(obj, ecell4_base.core.ReactionRule):
        for sp in itertools.chain(obj.reactants(), obj.products()):
            monomers = list_monomers(sp, monomers)
        return monomers
    elif isinstance(obj, ecell4_base.core.Species):
        for usp in obj.units():
            monomers = list_monomers(usp, monomers)
        return monomers
    elif not isinstance(obj, ecell4_base.core.UnitSpecies):
        raise TypeError(f"Invalid type [{repr(obj)}]")

    if obj.name() not in monomers:
        monomers[obj.name()] = {}
    for name, (state, bond) in obj.sites():
        if name not in monomers[obj.name()]:
            monomers[obj.name()][name] = []
        if state != '' and not state.startswith('_'):
            if state not in monomers[obj.name()][name]:  #XXX
                monomers[obj.name()][name].append(state)
    return monomers

def port_species(obj, model):
    if isinstance(obj, Iterable):
        return [port_species(obj_, model) for obj_ in obj]

    assert isinstance(obj, ecell4_base.core.Species)
    units = [port_unit_species(usp, model) for usp in obj.units()]
    assert len(units) != 0
    return functools.reduce(operator.mod, units)

def port_unit_species(obj, model):
    assert isinstance(obj, ecell4_base.core.UnitSpecies)
    opts = {}
    for name, (state, bond) in obj.sites():
        if state == '':
            if bond == '':
                opts[name] = None
            elif bond == '_':
                opts[name] = pysb.ANY
            elif bond.isdigit():
                opts[name] = int(bond)
            else:
                raise ValueError(f"Not supported {obj.serial()}")
        elif not state.startswith("_"):
            if bond == '':
                opts[name] = state
            elif bond == '_':
                opts[name] = (state, pysb.ANY)
            elif bond.isdigit():
                opts[name] = (state, int(bond))
            else:
                raise ValueError(f"Not supported {obj.serial()}")
    return model.monomers[obj.name()](**opts)

def convert(model, y0, name="MODEL", prefix="_", _export=False):
    assert isinstance(model, ecell4_base.core.NetfreeModel)  #XXX

    ret = PySBModel(name=name, _export=_export)

    for name, sites in list_monomers(model).items():
        ret.add_monomer(
            name,
            tuple(sites.keys()),
            {key: value for key, value in sites.items() if len(value) != 0})

    ignore = []
    for i, rr in enumerate(model.reaction_rules()):
        #XXX: Concat reversible reactions if exists. Avoid a problem in `pyvipr.pysb_viz.sp_rules_view`
        if i in ignore:
            continue
        kr = None
        rev = ecell4_base.core.ReactionRule(rr.products(), rr.reactants(), 0.0)
        for j, another in enumerate(model.reaction_rules()[i + 1: ]):
            if another == rr:
                raise RuntimeError("A reaction is duplicated [{}]".format(rr.as_string()))
            elif another == rev:
                assert (j + i + 1) not in ignore
                ignore.append(j + i + 1)
                if kr is not None:
                    raise RuntimeError("A reaction is duplicated [{}]".format(another.as_string()))
                kr = ret.add_parameter(f'{name}_kr', another.k())

        name = f'{prefix}RULE{i}'
        kf = ret.add_parameter(f'{name}_kf', rr.k())
        ret.add_rule(name, port_species(rr.reactants(), ret), port_species(rr.products(), ret), kf, kr)

    for i, (key, value) in enumerate(y0.items()):
        name = f'{prefix}INITIAL{i}'
        p = ret.add_parameter(name, value)
        ret.add_initial(port_species(ecell4_base.core.Species(key), ret), p)

    return ret
