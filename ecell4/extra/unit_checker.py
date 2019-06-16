from logging import getLogger
_log = getLogger(__name__)

import warnings
import operator
from functools import reduce

from ecell4_base.core import Quantity_Integer, Quantity_Real, Species

from ecell4.util import parseobj
from ecell4.util.decorator_base import just_parse
from ecell4.util.model_parser import Visitor, dispatch, RATELAW_RESERVED_FUNCTIONS

from ecell4.extra._unit import *


__all__ = ['check_units']

def is_structure(model, sp):
    sp = model.apply_species_attributes(sp)
    if not sp.has_attribute('structure'):
        return False
    val = sp.get_attribute('structure')
    assert isinstance(val, bool)
    return val

def assert_units(q, dim, key=''):
    q = q if isinstance(q, _Quantity) else getUnitRegistry().Quantity(q.magnitude, q.units)
    if not q.check(dim):
        raise RuntimeError('"{}" has a wrong dimensionality "{:s}". "{:s}" must be given.'.format(
            key, q.dimensionality, getUnitRegistry().get_dimensionality(dim)))

def get_dimension(model, sp):
    sp = model.apply_species_attributes(sp)
    if sp.has_attribute('dimension'):
        val = sp.get_attribute('dimension')
        assert isinstance(val, Quantity_Integer)
        return int(val.magnitude)
    elif sp.has_attribute('location'):
        loc = Species(sp.get_attribute('location'))
        return get_dimension(model, loc)
    return 3  # default

def get_units(model, sp):
    sp = model.apply_species_attributes(sp)
    if sp.has_attribute('units'):
        val = sp.get_attribute('units')
        assert isinstance(val, str)
        return val
    if is_structure(model, sp):
        ndim = get_dimension(model, sp)
        return "meter ** {:d}".format(ndim)
    return "item"  #XXX: Default units

class DimensionalityChecker(Visitor):

    OPERATORS = {
        parseobj.PosExp: operator.pos,
        parseobj.NegExp: operator.neg,
        parseobj.SubExp: lambda *args: operator.sub, # reduce(operator.sub, args[1: ], args[0]),
        parseobj.DivExp: operator.truediv, # lambda *args: reduce(operator.truediv, args[1: ], args[0]),
        parseobj.PowExp: operator.pow,
        parseobj.AddExp: lambda *args: reduce(operator.add, args[1: ], args[0]), # operator.add,
        parseobj.MulExp: lambda *args: reduce(operator.mul, args[1: ], args[0]), # operator.mul,
        # parseobj.InvExp: operator.inv,
        # parseobj.AndExp: operator.and_,
        # parseobj.GtExp: operator.gt,
        # parseobj.NeExp: operator.ne,
        # parseobj.EqExp: operator.eq,
        }

    def __init__(self, model, ureg):
        Visitor.__init__(self)
        self.__model = model
        self.__ureg = ureg

    def visit_const(self, obj):
        key = obj._elems[0].name
        if key == '_t':
            dim = self.__ureg.parse_units("second")
        elif key == '_v':
            dim = self.__ureg.parse_units("meter ** 3")
        else:
            dim = self.__ureg.parse_units("dimensionless")
        return dim

    def visit_species(self, obj):
        sp = Species(str(obj))
        dim = get_units(self.__model, sp)
        return self.__ureg.parse_units(dim)

    def visit_func(self, obj, *args):
        func = RATELAW_RESERVED_FUNCTIONS[obj._elems[0].name]
        val = func(*(1.0 * x for _, x in args))
        dim = val.to_base_units().u
        return dim

    def visit_expression(self, obj, *args):
        assert len(obj._elems) == len(args)
        for cls, op in self.OPERATORS.items():
            if isinstance(obj, cls):
                val = op(*[1.0 * x for x in args])
                dim = val.to_base_units().u
                return dim
        raise ValueError('Unknown dimensionality for the given object [{}]'.format(str(obj)))

    def visit_default(self, obj):
        if isinstance(obj, (Quantity_Integer, Quantity_Real)):
            dim = self.__ureg.parse_units(obj.units)
            return dim
        return obj

def check_units_ratelaw(model, rr):
    assert rr.has_descriptor()
    desc = rr.get_descriptor()
    name = desc.as_string()
    if name == '':
        warnings.warn('ReactionRuleDescriptor is not supported yet. "{}" was ignored.'.format(rr.as_string()))
        return

    obj = just_parse().eval(name, params={'Quantity_Integer': Quantity_Integer, 'Quantity_Real': Quantity_Real})
    ureg = getUnitRegistry()
    visitor = DimensionalityChecker(model, ureg)
    try:
        dim = dispatch(obj, visitor)
    except Exception as err:
        raise RuntimeError('An exception occurs while processing [{} ({})].'.format(name, rr.as_string()))

    dim = ureg.Quantity(1.0, dim).to_base_units().u
    assert_units(ureg.Quantity(1.0, dim), '[substance]/[time]', '{} ({})'.format(name, rr.as_string()))

def check_units(model):
    """
    Check units in the given model.

    Parameters
    ----------
    model : ecell4.core.Model
        A model to be checked.

    """
    dimensions = {
        'D': '[length]**2/[time]',
        'radius': '[length]',
        'dimension': '',  # means dimensionality of 'dimensionless'
    }

    for sp in model.species_attributes():
        for key, val in sp.list_attributes():
            if not isinstance(val, (Quantity_Integer, Quantity_Real)):
                continue
            elif val.units == '':
                continue
            elif key not in dimensions:
                continue
            assert_units(val, dimensions[key], '{} ({})'.format(key, sp.serial()))

    for rr in model.reaction_rules():
        if rr.has_descriptor():
            check_units_ratelaw(model, rr)
            continue

        val = rr.get_k()
        if val.units == '':
            continue

        reactants = rr.reactants()
        if len(reactants) == 0:
            products = rr.products()
            dim = set((get_dimension(model, sp) for sp in products))
            if len(dim) != 1:
                raise RuntimeError('All products must have the same dimension: "{}".'.format(rr.as_string()))
            basedim = int(dim[0])
            assert basedim in (1, 2, 3)
            assert_units(val, '[substance]/[time]/[length]**{:d}'.format(basedim), rr.as_string())
        else:
            dim = [get_dimension(model, sp) for sp in reactants]
            basedim = max(dim)
            structure = [is_structure(model, sp) for sp in reactants]
            units = '1/[time]/([substance]/[length]**{:d})**{:d}'.format(basedim, len(dim) - 1)
            if any(structure):
                tot = sum(dim_ for structure_, dim_ in zip(structure, dim) if structure_)
                units += '/([length]**{:d}/[substance]**{:d})'.format(tot, sum(structure))
            assert_units(val, units, rr.as_string())
