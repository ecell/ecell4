import numbers
from .session import load_world, Session

from .decorator import get_model, reset_model
from . import viz
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
        ``ODEWorld``, ``GillespieWorld``, ``SpatiocyteWorld`` and ``SGFRDWorld``.

    """
    import ecell4_base

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
    elif vinfo.startswith("ecell4-sgfrd"):
        return ecell4_base.sgfrd.World(filename)
    elif vinfo == "":
        raise RuntimeError("No version information was found in [{0}]".format(filename))
    raise RuntimeError("Unknown version information [{0}]".format(vinfo))

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
    elif solver == 'sgfrd':
        return ecell4_base.sgfrd.Factory(*args)
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
        t, y0=None, volume=1.0, model=None, solver='ode', ndiv=None,
        species_list=None, structures=None, observers=(), rndseed=None):
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
        'bd', 'egfrd', and 'sgfrd'. Default is 'ode'.
        When tuple is given, the first value must be str as explained above.
        All the rest is used as arguments for the corresponding factory class.
        Keyword 's' is a shortcut for specifying 'solver'.
    ndiv : int, optional
        A number of time points. If t is an array, ignored. If None, log all.
    species_list : list of str, optional
        A list of names of Species observed. If None, log all.
        Default is None.
    structures : list or dict, optional
        A dictionary which gives pairs of a name and shape of structures.
        Not fully supported yet.
    observers : Observer or list, optional
        A list of extra observer references.
    rndseed : int, optional
        A random seed for a simulation.
        This argument will be ignored when 'solver' is given NOT as a string.

    Returns
    -------
    result : Result

    """
    session = Session(model=model, y0=y0, structures=structures, volume=volume)
    ret = session.run(
        t, solver=solver, rndseed=rndseed, ndiv=ndiv, species_list=species_list, observers=observers)
    return ret

def ensemble_simulations(
        t, y0=None, volume=1.0, model=None, solver='ode', ndiv=None,
        species_list=None, structures=None, observers=(), rndseed=None,
        repeat=1, method=None,
        **kwargs):
    """
    Run simulations multiple times and return its ensemble.
    Arguments are almost same with ``ecell4.util.simulation.run_simulation``.
    `observers` and `progressbar` is not available here.

    Parameters
    ----------
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
    value : ResultList

    See Also
    --------
    ecell4.util.simulation.run_simulation
    ecell4.extra.ensemble.run_serial
    ecell4.extra.ensemble.run_sge
    ecell4.extra.ensemble.run_slurm
    ecell4.extra.ensemble.run_multiprocessing
    ecell4.extra.ensemble.run_azure

    """
    session = Session(model=model, y0=y0, structures=structures, volume=volume)
    ret = session.ensemble(
        t, solver=solver, rndseed=rndseed, ndiv=ndiv, species_list=species_list, observers=observers,
        repeat=repeat, method=method, **kwargs)
    return ret

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
    from ecell4_base.core import NumberObserver, FixedIntervalNumberObserver, TimingNumberObserver

    if t is None:
        return NumberObserver(targets)
    elif isinstance(t, numbers.Number):
        return FixedIntervalNumberObserver(t, targets)
    elif hasattr(t, '__iter__'):
        if targets is not None:
            return TimingNumberObserver(t, targets)
        else:
            return TimingNumberObserver(t)
    else:
        raise TypeError("An invalid type was given. Either number or iterable is required.")
