import numbers
from .session import load_world, Session


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
        'bd' and 'egfrd'. Default is 'ode'.
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
