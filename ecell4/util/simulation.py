import collections
import numbers

from .decorator import get_model, reset_model
from . import viz
from ..extra import unit
from ..extra.ensemble import ensemble_simulations


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
    import ecell4_base

    vinfo = ecell4_base.core.load_version_information(filename)
    if vinfo.startswith("ecell4-bd"):
        return ecell4.bd.World(filename)
    elif vinfo.startswith("ecell4-egfrd"):
        return ecell4.egfrd.World(filename)
    elif vinfo.startswith("ecell4-meso"):
        return ecell4.meso.World(filename)
    elif vinfo.startswith("ecell4-ode"):
        return ecell4.ode.World(filename)
    elif vinfo.startswith("ecell4-gillespie"):
        return ecell4.gillespie.World(filename)
    elif vinfo.startswith("ecell4-spatiocyte"):
        return ecell4.spatiocyte.World(filename)
    elif vinfo == "":
        raise RuntimeError("No version information was found in [{0}]".format(filename))
    raise RuntimeError("Unkown version information [{0}]".format(vinfo))

def get_factory(solver, *args):
    import ecell4_base

    if solver == 'ode':
        return ecell4.ode.Factory(*args)
    elif solver == 'gillespie':
        return ecell4.gillespie.Factory(*args)
    elif solver == 'spatiocyte':
        return ecell4.spatiocyte.Factory(*args)
    elif solver == 'meso':
        return ecell4.meso.Factory(*args)
    elif solver == 'bd':
        return ecell4.bd.Factory(*args)
    elif solver == 'egfrd':
        return ecell4.egfrd.Factory(*args)
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

def run_simulation(
        t, y0=None, volume=1.0, model=None, solver='ode',
        is_netfree=False, species_list=None, without_reset=False,
        return_type='matplotlib', opt_args=(), opt_kwargs=None,
        structures=None, observers=(), progressbar=0, rndseed=None,
        factory=None, ## deprecated
        **kwargs):
    """Run a simulation with the given model and plot the result on IPython
    notebook with matplotlib.

    Parameters
    ----------
    t : array or Real
        A sequence of time points for which to solve for 'm'.
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.
        Keyword 'v' is a shortcut for specifying 'volume'.
    model : Model, optional
        Keyword 'm' is a shortcut for specifying 'model'.
    solver : str, tuple or Factory, optional
        Solver type. Choose one from 'ode', 'gillespie', 'spatiocyte', 'meso',
        'bd' and 'egfrd'. Default is 'ode'.
        When tuple is given, the first value must be str as explained above.
        All the rest is used as arguments for the corresponding factory class.
        Keyword 's' is a shortcut for specifying 'solver'.
    species_list : list of str, optional
        A list of names of Species observed. If None, log all.
        Default is None.
    return_type : str, optional
        Choose a type of return value from 'array', 'observer',
        'matplotlib', 'nyaplot', 'world', 'dataframe', 'none' or None.
        If None or 'none', return and plot nothing. Default is 'matplotlib'.
        'dataframe' requires numpy and pandas libraries.
        Keyword 'r' is a shortcut for specifying 'return_type'.
    opt_args: list, tuple or dict, optional
        Arguments for plotting. If return_type suggests no plotting, just ignored.
    opt_kwargs: dict, optional
        Arguments for plotting. If return_type suggests no plotting or
        opt_args is a list or tuple, just ignored.
        i.e.) viz.plot_number_observer(obs, *opt_args, **opt_kwargs)
    is_netfree: bool, optional
        Whether the model is netfree or not. When a model is given as an
        argument, just ignored. Default is False.
    structures : list or dict, optional
        A dictionary which gives pairs of a name and shape of structures.
        Not fully supported yet.
    observers : Observer or list, optional
        A list of extra observer references.
    progressbar : float, optional
        A timeout for a progress bar in seconds.
        When the value is not more than 0, show nothing.
        Default is 0.
    rndseed : int, optional
        A random seed for a simulation.
        This argument will be ignored when 'solver' is given NOT as a string.

    Returns
    -------
    value : list, TimingNumberObserver, World or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return an observer.
        When ``return_type`` is 'world', return the last state of ``World``.
        Return nothing if else.

    """
    y0 = y0 or {}
    opt_kwargs = opt_kwargs or {}
    structures = structures or {}

    for key, value in kwargs.items():
        if key == 'r':
            return_type = value
        elif key == 'v':
            volume = value
        elif key == 's':
            solver = value
        elif key == 'm':
            model = value
        else:
            raise ValueError(
                "An unknown keyword argument was given [{}={}]".format(key, value))

    import ecell4_base

    if unit.HAS_PINT:
        if isinstance(t, unit._Quantity):
            if unit.STRICT and not unit.check_dimensionality(t, '[time]'):
                raise ValueError("Cannot convert [t] from '{}' ({}) to '[time]'".format(t.dimensionality, t.u))
            t = t.to_base_units().magnitude

        if isinstance(volume, unit._Quantity):
            if unit.STRICT:
                if isinstance(volume.magnitude, ecell4_base.core.Real3) and not unit.check_dimensionality(volume, '[length]'):
                    raise ValueError("Cannot convert [volume] from '{}' ({}) to '[length]'".format(
                        volume.dimensionality, volume.u))
                elif not unit.check_dimensionality(volume, '[volume]'):
                    raise ValueError("Cannot convert [volume] from '{}' ({}) to '[volume]'".format(
                        volume.dimensionality, volume.u))
            volume = volume.to_base_units().magnitude

        if not isinstance(solver, str) and isinstance(solver, collections.Iterable):
            solver = [
                value.to_base_units().magnitude if isinstance(value, unit._Quantity) else value
                for value in solver]

    if factory is not None:
        # f = factory  #XXX: will be deprecated in the future. just use solver
        raise ValueError(
            "Argument 'factory' is no longer available. Use 'solver' instead.")
    elif isinstance(solver, str):
        f = get_factory(solver)
    elif isinstance(solver, collections.Iterable):
        f = get_factory(*solver)
    else:
        f = solver

    if rndseed is not None:
        f = f.rng(ecell4_base.core.GSLRandomNumberGenerator(rndseed))

    if model is None:
        model = ecell4.util.decorator.get_model(is_netfree, without_reset)

    w = f.world(volume)
    edge_lengths = w.edge_lengths()

    if unit.HAS_PINT:
        y0 = y0.copy()
        for key, value in y0.items():
            if isinstance(value, unit._Quantity):
                if not unit.STRICT:
                    y0[key] = value.to_base_units().magnitude
                elif unit.check_dimensionality(value, '[substance]'):
                    y0[key] = value.to_base_units().magnitude
                elif unit.check_dimensionality(value, '[concentration]'):
                    volume = w.volume() if not isinstance(w, ecell4.spatiocyte.SpatiocyteWorld) else w.actual_volume()
                    y0[key] = value.to_base_units().magnitude * volume
                else:
                    raise ValueError(
                        "Cannot convert a quantity for [{}] from '{}' ({}) to '[substance]'".format(
                            key, value.dimensionality, value.u))

    if not isinstance(w, ecell4.ode.ODEWorld):
        w.bind_to(model)

    for (name, shape) in (structures.items() if isinstance(structures, dict) else structures):
        if isinstance(shape, str):
            w.add_structure(ecell4_base.core.Species(name), get_shape(shape))
        elif isinstance(shape, collections.Iterable):
            w.add_structure(ecell4_base.core.Species(name), get_shape(*shape))
        else:
            w.add_structure(ecell4_base.core.Species(name), shape)

    if isinstance(w, ecell4.ode.ODEWorld):
        # w.bind_to(model)  # stop binding for ode
        for serial, n in y0.items():
            w.set_value(ecell4_base.core.Species(serial), n)
    else:
        # w.bind_to(model)
        for serial, n in y0.items():
            w.add_molecules(ecell4_base.core.Species(serial), n)

    if not isinstance(t, collections.Iterable):
        t = [float(t) * i / 100 for i in range(101)]

    obs = ecell4_base.core.TimingNumberObserver(t, species_list)
    sim = f.simulator(w, model)
    # sim = f.simulator(w)

    if not isinstance(observers, collections.Iterable):
        observers = (observers, )
    if return_type not in ('world', 'none', None):
        observers = (obs, ) + tuple(observers)

    if progressbar > 0:
        from .progressbar import progressbar as pb
        pb(sim, timeout=progressbar, flush=True).run(t[-1], observers)
    else:
        sim.run(t[-1], observers)

    if return_type in ('matplotlib', 'm'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.viz.plot_number_observer(obs, *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.viz.plot_number_observer(obs, **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('nyaplot', 'n'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.viz.plot_number_observer_with_nya(obs, *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.viz.plot_number_observer_with_nya(obs, **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('observer', 'o'):
        return obs
    elif return_type in ('array', 'a'):
        return obs.data()
    elif return_type in ('dataframe', 'd'):
        import pandas
        import numpy
        data = numpy.array(obs.data()).T
        return pandas.concat([
            pandas.DataFrame(dict(Time=data[0], Value=data[i + 1],
                                  Species=sp.serial(), **opt_kwargs))
            for i, sp in enumerate(obs.targets())])
    elif return_type in ('world', 'w'):
        return sim.world()
    elif return_type is None or return_type in ('none', ):
        return
    else:
        raise ValueError(
            'An invald value for "return_type" was given [{}].'.format(str(return_type))
            + 'Use "none" if you need nothing to be returned.')

def number_observer(t=None, targets=None):
    """
    Return a number observer. If t is None, return NumberObserver. If t is a number,
    return FixedIntervalNumberObserver. If t is an iterable (a list of numbers), return
    TimingNumberObserver.

    Parameters
    ----------
    t : float, list or tuple, optional. default None
        A timing of the observation. See above.
    targets : list or tuple, optional. default None
        A list of strings suggesting Species observed.

    Returns
    -------
    obs : NumberObserver, FixedIntervalNumberObserver or TimingNumberObserver
    """
    from ecell4 import NumberObserver, FixedIntervalNumberObserver, TimingNumberObserver

    if t is None:
        return NumberObserver(targets)
    elif isinstance(t, numbers.Number):
        return FixedIntervalNumberObserver(t, targets)
    elif hasattr(t, '__iter__'):
        return TimingNumberObserver(t, targets)
    else:
        raise TypeError("An invalid type was given. Either number or iterable is required.")
