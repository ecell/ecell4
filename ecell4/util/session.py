import copy
import collections
import numbers
import ecell4_base

from . import viz
from .decorator import get_model
from ..extra import unit


def get_factory(solver, *args):
    import ecell4_base

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

    import ecell4_base
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

        #XXX: y0 may accept [concentration].
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

    def ensemble(
        self, t, solver='ode', rndseed=None, ndiv=None, species_list=None,
        repeat=1, nproc=None, method=None, **kwargs):
        """
        Run simulations multiple times and return its ensemble.
        Arguments are almost same with ``ecell4.util.simulation.run_simulation``.
        `observers` and `progressbar` is not available here.

        Parameters
        ----------
        repeat : int, optional
            A number of runs. Default is 1.
        nproc : int, optional
            A number of processors. Ignored when method='serial'.
            Default is None.
        method : str, optional
            The way for running multiple jobs.
            Choose one from 'serial', 'multiprocessing', 'sge', 'slurm', 'azure'.
            Default is None, which works as 'serial'.
        **kwargs : dict, optional
            Optional keyword arugments are passed through to `run_serial`,
            `run_sge`, or `run_multiprocessing`.
            See each function for more details.

        See Also
        --------
        ecell4.extra.ensemble.run_serial
        ecell4.extra.ensemble.run_sge
        ecell4.extra.ensemble.run_slurm
        ecell4.extra.ensemble.run_multiprocessing
        ecell4.extra.ensemble.run_azure

        """
        from ..extra.ensemble import genseeds
        # if species_list is None:
        #     species_list = list_species(model, y0.keys())

        if rndseed is None:
            myseed = genseeds(repeat)
        elif (not isinstance(rndseed, bytes) or len(rndseed) != repeat * 4 * 2):
            raise ValueError(
                "A wrong seed for the random number generation was given. Use 'genseeds'.")

        jobs = [{'session': self, 't': t, 'solver': solver, 'myseed': myseed, 'ndiv': ndiv, 'species_list': species_list}]

        from ..extra.ensemble import run_serial, run_sge, run_slurm, run_multiprocessing, run_azure
        if method is None or method.lower() == "serial":
            results = run_serial(singlerun, jobs, n=repeat, **kwargs)
        elif method.lower() == "sge":
            results = run_sge(singlerun, jobs, n=repeat, nproc=nproc, **kwargs)
        elif method.lower() == "slurm":
            results = run_slurm(singlerun, jobs, n=repeat, nproc=nproc, **kwargs)
        elif method.lower() == "multiprocessing":
            results = run_multiprocessing(singlerun, jobs, n=repeat, nproc=nproc, **kwargs)
        elif method.lower() == "azure":
            results = run_azure(singlerun, jobs, n=repeat, nproc=nproc, **kwargs)
        else:
            raise ValueError(
                'Argument "method" must be one of "serial", "multiprocessing", "slurm" and "sge".')

        assert len(results) == len(jobs) == 1
        assert len(results[0]) == repeat
        return results
