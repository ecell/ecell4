import warnings
import copy
import collections.abc
import numbers
import tempfile
import os

import ecell4_base

from .decorator import get_model
from ..extra import unit


def load_world(filename):
    """
    Load a world from the given HDF5 filename.
    The return type is determined by ``ecell4_base.core.load_version_information``.

    Parameters
    ----------
    filename : str
        A HDF5 filename.

    Returns
    -------
    w : World
        Return one from ``BDWorld``, ``EGFRDWorld``, ``MesoscopicWorld``,
        ``ODEWorld``, ``GillespieWorld`` and ``SpatiocyteWorld``.

    """
    vinfo = ecell4_base.core.load_version_information(filename)
    if vinfo.startswith("ecell4-bd"):
        return ecell4_base.bd.World(filename)
    elif vinfo.startswith("ecell4-egfrd"):
        return ecell4_base.egfrd.World(filename)
    elif vinfo.startswith("ecell4-meso"):
        return ecell4_base.meso.World(filename)
    elif vinfo.startswith("ecell4-ode"):
        return ecell4_base.ode.World(filename)
    elif vinfo.startswith("ecell4-gillespie"):
        return ecell4_base.gillespie.World(filename)
    elif vinfo.startswith("ecell4-spatiocyte"):
        return ecell4_base.spatiocyte.World(filename)
    elif vinfo == "":
        raise RuntimeError("No version information was found in [{0}]".format(filename))
    raise RuntimeError("Unknown version information [{0}]".format(vinfo))

def get_factory(solver, *args):
    if solver == 'ode':
        return ecell4_base.ode.Factory(*args)
    elif solver == 'gillespie':
        return ecell4_base.gillespie.Factory(*args)
    elif solver == 'spatiocyte':
        return ecell4_base.spatiocyte.Factory(*args)
    elif solver == 'meso':
        return ecell4_base.meso.Factory(*args)
    elif solver == 'bd':
        return ecell4_base.bd.Factory(*args)
    elif solver == 'egfrd':
        return ecell4_base.egfrd.Factory(*args)
    else:
        raise ValueError(
            'unknown solver name was given: ' + repr(solver)
            + '. use ode, gillespie, spatiocyte, meso, bd or egfrd')

def get_shape(shape, *args):
    if not isinstance(shape, str):
        raise ValueError("Invalid shape was given [{}]. This must be 'str'".format(repr(shape)))

    shape = shape.lower()
    shape_map = {
        'aabb': ecell4_base.core.AABB,
        # 'affinetransformation': ecell4_base.core.AffineTransformation,
        'cylinder': ecell4_base.core.Cylinder,
        'cylindricalsurface': ecell4_base.core.CylindricalSurface,
        'meshsurface': ecell4_base.core.MeshSurface,
        'planarsurface': ecell4_base.core.PlanarSurface,
        'rod': ecell4_base.core.Rod,
        'rodsurface': ecell4_base.core.RodSurface,
        'sphere': ecell4_base.core.Sphere,
        'sphericalsurface': ecell4_base.core.SphericalSurface,
        # 'complement': ecell4_base.core.Complement,
        # 'union': ecell4_base.core.Union,
        }
    if shape in shape_map:
        args = [
            value.to_base_units().magnitude if isinstance(value, unit._Quantity) else value
            for value in args]
        return shape_map[shape](*args)
    else:
        raise ValueError(
            'unknown shape type was given: ' + repr(shape)
            + '. use {}'.format(', '.join(sorted(shape_map.keys()))))

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

def singlerun(job, job_id, task_id):
    import ecell4.extra.ensemble
    session = job.pop('session')
    rndseed = ecell4.extra.ensemble.getseed(job.pop('myseed'), task_id)
    job['rndseed'] = rndseed
    return session.run(**job)

class Result(object):

    def __init__(self, world, observers):
        self.world = world
        self.observers = observers

    @property
    def observer(self):
        return self.observers[0]

    # @property
    # def y0(self):
    #     return [(sp.serial(), self.world.get_value_exact(sp)) for sp in self.world.list_species()]

    # @property
    # def volume(self):
    #     return self.world.edge_lengths()

    def data(self):
        return self.observer.data()

    def targets(self):
        return self.observer.targets()

    def species_list(self):
        return [sp.serial() for sp in self.observer.targets()]

    def plot(self, *args, **kwargs):
        from ..plotting import plot_number_observer
        plot_number_observer(self.observer, *args, **kwargs)

    def as_array(self, *args, **kwargs):
        """Require numpy"""
        import numpy
        return numpy.array(self.observer.data(), *args, **kwargs)

    def as_df(self):
        """See as_dataframe."""
        return self.as_dataframe()

    def as_dataframe(self):
        """Require pandas"""
        import pandas
        return pandas.DataFrame(data=self.data(), columns=['t'] + self.species_list())

    def __getstate__(self):
        fd, tmpfile = tempfile.mkstemp()
        try:
            self.world.save(tmpfile)
            with open(tmpfile, 'rb') as fp:
                world = fp.read()
        finally:
            os.remove(tmpfile)
        return (world, self.observers)

    def __setstate__(self, state):
        if len(state) != 2:
            raise ValueError("Invalid state!")
        fd, tmpfile = tempfile.mkstemp()
        try:
            with open(tmpfile, 'wb') as fp:
                fp.write(state[0])
            world = load_world(tmpfile)
        finally:
            os.remove(tmpfile)
        self.world = world
        self.observers = state[1]

    def _ipython_display_(self):
        """
        Displays the object as a side effect.
        https://ipython.readthedocs.io/en/stable/config/integrating.html

        """
        self.plot()

class ResultList(object):

    def __init__(self, items):
        for item in items:
            if not isinstance(item, Result):
                raise ValueError("The item must be a Result object [{}].".format(repr(item)))
        self.__items = items
        self.__initialize()

    def __initialize(self):
        self.__data = []
        self.__err = []
        self.__targets = []

        if len(self.__items) == 0:
            warnings.warn("No item was given.")
            return
        elif len(self.__items) == 1:
            warnings.warn("Only one item was given.")
            self.__data = self.items[0].data()
            self.__err = []
            self.__targets = self.items[0].targets()
            return

        import numpy
        t = self.__items[0].as_array(numpy.float64).T[0]
        if len(t) == 0:
            warnings.warn("No data was given.")
            return

        data = []
        for item in self.__items:
            data_ = item.as_array(numpy.float64).T[1: ]
            if data_.shape[1] != len(t):
                warnings.warn("Invalid state. The length of data varies [{} != {}].".format(len(t), data_.shape[1]))
                return
            data.append(data_)

        self.__targets = self.__items[0].targets()
        self.__data = numpy.vstack([t, numpy.average(data, axis=0)]).T
        self.__err = numpy.vstack([t, numpy.std(data, axis=0)]).T

    def data(self):
        return self.__data

    def error(self):
        return self.__err

    def targets(self):
        return self.__targets

    def species_list(self):
        return [sp.serial() for sp in self.__targets]

    def plot(self, *args, **kwargs):
        from ..plotting import plot_number_observer
        plot_number_observer(self, *args, **kwargs)

    # def as_array(self):
    #     """Require numpy"""
    #     import numpy
    #     return numpy.array(self.observer.data())

    # def as_df(self):
    #     """See as_dataframe."""
    #     return self.as_dataframe()

    # def as_dataframe(self):
    #     """Require pandas"""
    #     import pandas
    #     return pandas.DataFrame(data=self.data(), columns=['t'] + self.species_list())

    def __getstate__(self):
        return (self.__items, )

    def __setstate__(self, state):
        if len(state) != 1:
            raise ValueError("Invalid state!")
        self.__items = state[0]
        self.__initialize()

    def _ipython_display_(self):
        """
        Displays the object as a side effect.
        https://ipython.readthedocs.io/en/stable/config/integrating.html

        """
        self.plot()

    def __len__(self):
        return len(self.__items)

    def __getitem__(self, key):
        return self.__items.__getitem__(key)

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
        import numpy
        t = get_value(t, ((numbers.Real, '[time]'), (collections.abc.Iterable, '[time]')))

        if isinstance(solver, str):
            f = get_factory(solver)
        elif isinstance(solver, collections.abc.Iterable):
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
            elif isinstance(shape, collections.abc.Iterable):
                w.add_structure(ecell4_base.core.Species(name), get_shape(*shape))
            else:
                w.add_structure(ecell4_base.core.Species(name), shape)  # shape is a Shape

        #XXX: y0 may accept [concentration].
        if isinstance(w, ecell4_base.ode.ODEWorld):
            for serial, n in self.y0.items():
                w.set_value(ecell4_base.core.format_species(ecell4_base.core.Species(serial)), n)
        else:
            for serial, n in self.y0.items():
                w.add_molecules(ecell4_base.core.format_species(ecell4_base.core.Species(serial)), int(n))

        if isinstance(w, ecell4_base.ode.ODEWorld):
            ndiv = ndiv or 100

        if isinstance(t, collections.abc.Iterable):
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

        if not isinstance(observers, collections.abc.Iterable):
            observers = (observers, )
        observers = (obs, ) + tuple(observers)

        sim = f.simulator(w, self.model)
        sim.run(upto, observers)

        return Result(w, observers)

    def ensemble(
        self, t, solver='ode', rndseed=None, ndiv=None, species_list=None, observers=(),
        repeat=1, method=None, **kwargs):
        """
        Run simulations multiple times and return its ensemble.
        Arguments are almost same with ``ecell4.util.simulation.run_simulation``.
        `observers` and `progressbar` is not available here.

        Parameters
        ----------
        t : array or Real
            A sequence of time points for which to solve for 'm'.
        solver : str, tuple or Factory, optional
            Solver type. Choose one from 'ode', 'gillespie', 'spatiocyte', 'meso',
            'bd' and 'egfrd'. Default is 'ode'.
            When tuple is given, the first value must be str as explained above.
            All the rest is used as arguments for the corresponding factory class.
        ndiv : int, optional
            A number of time points. If t is an array, ignored. If None, log all.
        species_list : list of str, optional
            A list of names of Species observed. If None, log all.
            Default is None.
        observers : Observer or list, optional
            A list of extra observer references.
        repeat : int, optional
            A number of runs. Default is 1.
        method : str, optional
            The way for running multiple jobs.
            Choose one from 'serial', 'multiprocessing', 'sge', 'slurm', 'azure'.
            Default is None, which works as 'serial'.
        **kwargs : dict, optional
            Optional keyword arugments are passed through to `run_serial`,
            `run_sge`, or `run_multiprocessing`.
            See each function for more details.

        Returns
        -------
        results : list
            A list of Result objects. The list contains `n` results.

        See Also
        --------
        ecell4.extra.ensemble.run_ensemble

        """
        from ..extra.ensemble import genseeds
        # if species_list is None:
        #     species_list = list_species(model, y0.keys())

        if rndseed is None:
            myseed = genseeds(repeat)
        elif (not isinstance(rndseed, bytes) or len(rndseed) != repeat * 4 * 2):
            raise ValueError(
                "A wrong seed for the random number generation was given. Use 'genseeds'.")

        jobs = [{
            'session': self, 't': t, 'solver': solver, 'myseed': myseed, 'ndiv': ndiv,
            'species_list': species_list, 'observers': observers,
            }]

        from ..extra.ensemble import run_ensemble
        results = run_ensemble(singlerun, jobs, repeat=repeat, method=method, **kwargs)

        assert len(results) == len(jobs) == 1
        assert len(results[0]) == repeat
        return ResultList(results[0])
