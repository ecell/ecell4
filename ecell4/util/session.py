import copy
import collections
import numbers
import ecell4_base

from . import viz
from .decorator import get_model
from .simulation import get_factory, get_shape, number_observer

from ..extra import unit


def get_value(obj, type_and_dim=None):
    if unit.HAS_PINT and isinstance(obj, unit._Quantity):
        for typeinfo, dim in type_and_dim:
            if isinstance(obj.magnitude, typeinfo) and obj.check(dim):
                return obj.to_base_units().magnitude
        else:
            raise TypeError('The given value [{}] has wrong type or dimension.'.format(repr(obj)))
    elif isinstance(obj, (ecell4_base.core.Quantity_Real, ecell4_base.core.Quantity_Integer)):
        return get_value(obj.magnitude, type_and_dim)

    if type_and_dim is not None:
        for typeinfo, dim in type_and_dim:
            if isinstance(obj, typeinfo):
                break
        else:
            raise TypeError('The given value [{}] has wrong type or dimension.'.format(repr(obj)))
    return obj

class Result(object):

    def __init__(self, world, observers):
        self.world = world
        self.observers = observers

    @property
    def observer(self):
        return self.observers[0]

    def plot(self, *args, **kwargs):
        viz.plot_number_observer(self.observer, *args, **kwargs)

    def data(self):
        return self.observer.data()

    def species_list(self):
        return [sp.serial() for sp in self.observer.targets()]

    def as_array(self):
        """Require numpy"""
        import numpy
        return numpy.array(self.observer.data())

    def as_df(self):
        return self.as_dataframe()

    def as_dataframe(self):
        """Require pandas"""
        import pandas
        return pandas.DataFrame(data=self.data(), columns=['t'] + self.species_list())

class Session(object):

    def __init__(self, model=None, y0=None, structures=None, volume=1.0):
        """
        Constructor

        Parameters
        ----------
        model : Model, optional
        y0 : dict
            Initial condition.
        structures : dict, optional
            A dictionary which gives pairs of a name and shape of structures.
            Not fully supported yet.
        volume : Real or Real3, optional
            A size of the simulation volume.
            1.0 as default.

        """
        self.model = model or get_model()

        self.y0 = {}
        if y0 is not None:
            for key, value in y0.items():
                self.y0[key] = get_value(value, ((numbers.Real, '[substance]'), ))

        self.structures = structures.copy() if structures is not None else {}
        self.volume = get_value(volume, ((numbers.Real, '[volume]'), (ecell4_base.core.Real3, '[length]')))

    def run(self, t, solver='ode', rndseed=None, ndiv=None, species_list=None, observers=()):
        """Run a simulation with the given model and return the result

        Parameters
        ----------
        t : array or Real
            A sequence of time points for which to solve for 'm'.
        solver : str, tuple or Factory, optional
            Solver type. Choose one from 'ode', 'gillespie', 'spatiocyte', 'meso',
            'bd' and 'egfrd'. Default is 'ode'.
            When tuple is given, the first value must be str as explained above.
            All the rest is used as arguments for the corresponding factory class.
        rndseed : int, optional
            A random seed for a simulation.
        ndiv : int, optional
            A number of time points. If t is an array, ignored. If None, log all.
        species_list : list of str, optional
            A list of names of Species observed. If None, log all.
            Default is None.
        observers : Observer or list, optional
            A list of extra observer references.

        Returns
        -------
        ret : Result
            The result object

        """

        t = get_value(t, ((numbers.Real, '[time]'), ))

        if isinstance(solver, str):
            f = get_factory(solver)
        elif isinstance(solver, collections.Iterable):
            f = get_factory(*solver)
        else:
            f = solver  # solver is a Factory

        if rndseed is not None:
            f = f.rng(ecell4_base.core.GSLRandomNumberGenerator(rndseed))

        w = f.world(self.volume)

        if not isinstance(w, ecell4_base.ode.ODEWorld):
            w.bind_to(self.model)

        for (name, shape) in self.structures.items():
            if isinstance(shape, str):
                w.add_structure(ecell4_base.core.Species(name), get_shape(shape))
            elif isinstance(shape, collections.Iterable):
                w.add_structure(ecell4_base.core.Species(name), get_shape(*shape))
            else:
                w.add_structure(ecell4_base.core.Species(name), shape)  # shape is a Shape

        if isinstance(w, ecell4_base.ode.ODEWorld):
            for serial, n in self.y0.items():
                w.set_value(ecell4_base.core.Species(serial), n)
        else:
            for serial, n in self.y0.items():
                w.add_molecules(ecell4_base.core.Species(serial), n)

        if isinstance(w, ecell4_base.ode.ODEWorld):
            ndiv = ndiv or 100

        if isinstance(t, collections.Iterable):
            upto = t[-1]
            if species_list is None:
                obs = ecell4_base.core.TimingNumberObserver(t)
            else:
                obs = ecell4_base.core.TimingNumberObserver(t, species_list)
        elif ndiv is not None:
            upto = t
            assert isinstance(ndiv, numbers.Integral) and int(ndiv) > 0
            t = [float(t) * i / int(ndiv) for i in range(int(ndiv) + 1)]
            if species_list is None:
                obs = ecell4_base.core.TimingNumberObserver(t)
            else:
                obs = ecell4_base.core.TimingNumberObserver(t, species_list)
        else:
            upto = t
            if species_list is None:
                obs = ecell4_base.core.NumberObserver()
            else:
                obs = ecell4_base.core.NumberObserver(species_list)

        if not isinstance(observers, collections.Iterable):
            observers = (observers, )
        observers = (obs, ) + tuple(observers)

        sim = f.simulator(w, self.model)
        sim.run(upto, observers)

        return Result(w, observers)
