import collections
import numbers

from .decorator import get_model, reset_model
from . import viz
from ..extra import unit
from ..extra.ensemble import ensemble_simulations
from .session import load_world


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

    from .session import Session
    session = Session(model=model, y0=y0, structures=structures, volume=volume)
    ret = session.run(t, solver=solver, rndseed=rndseed, ndiv=None, species_list=species_list, observers=observers)

    if return_type in ('matplotlib', 'm'):
        if isinstance(opt_args, (list, tuple)):
            ret.plot(*opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ret.plot(*opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('nyaplot', 'n'):
        if isinstance(opt_args, (list, tuple)):
            viz.plot_number_observer_with_nya(ret.observer, *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            viz.plot_number_observer_with_nya(ret.obs, **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('observer', 'o'):
        return ret.observer
    elif return_type in ('array', 'a'):
        return ret.data()
    elif return_type in ('dataframe', 'd'):
        return ret.as_df()
    elif return_type in ('world', 'w'):
        return ret.world
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
