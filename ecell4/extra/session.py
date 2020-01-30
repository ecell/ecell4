import numbers
import collections

import ecell4_base.core

from ecell4.util.decorator import get_model
from ecell4.util.simulation import get_factory
from ecell4.util import viz


class Result(object):

    def __init__(self, world, observers):
        self.world = world
        self.observers = observers

    @property
    def observer(self):
        return self.observers[0]

    @property
    def columns(self):
        return ['t'] + [sp.serial() for sp in self.observer.targets()]

    @property
    def values(self):
        import numpy
        return numpy.array(self.observer.data())

    def as_array(self):
        return self.values

    def as_dataframe(self):
        import pandas
        return pandas.DataFrame(self.values, columns=self.columns)

    def plot(self, *args, **kwargs):
        viz.plot_number_observer(self.observer, *args, **kwargs)

    @property
    def y0(self):
        return [(sp.serial(), self.world.get_value_exact(sp)) for sp in self.world.list_species()]

    @property
    def volume(self):
        return self.world.edge_lengths()

def run_tasks(tasks):
    return [run(**tasks) for task in tasks]

def run(t, y0=None, volume=1.0, model=None, solver='ode', rndseed=None):
    """

    Parameters
    ----------
    t : Real or iterable object of Real
        A sequence of time points.
    y0 : iterable object, optional
        Initial condition (the default is []).
    volume : Real or Real3, optional
        A size of the simulation volume (the default is 1).
    model : Model, optional
        The model to be solved. If `model` is `None`, call `ecell4.util.get_model`.
    solver : str, optional
        Solver kind. Choose one from 'ode', 'gillespie', 'spatiocyte', 'meso',
        'bd' and 'egfrd' (the default is 'ode').
        When tuple is given, the first value must be str as explained above.
        All the rest is used as arguments for the corresponding factory class.
    rndseed : int, optional
        A random seed.

    """
    y0 = y0 or []
    model = model or get_model()

    assert isinstance(y0, collections.Iterable)
    assert isinstance(volume, (numbers.Real, ecell4_base.core.Real3))
    assert isinstance(model, ecell4_base.core.Model)

    if isinstance(solver, str):
        f = get_factory(solver)
    elif isinstance(solver, collections.Iterable):
        f = get_factory(*solver)
    else:
        raise TypeError("'solver' must be str or iterable.")

    if isinstance(rndseed, numbers.Integral):
        f = f.rng(ecell4_base.core.GSLRandomNumberGenerator(rndseed))
    else:
        assert rndseed is None

    w = f.world(volume)
    w.bind_to(model)

    for key, val in y0:
        sp = ecell4_base.core.Species(key)
        if isinstance(val, numbers.Real):
            w.set_value(sp, val)
        else:
            w.add_structure(sp, val)

    if not isinstance(t, collections.Iterable):
        assert isinstance(t, numbers.Real)
        t = [float(t) * i / 100 for i in range(101)]

    obs = (ecell4_base.core.TimingNumberObserver(t), )

    sim = f.simulator(w, model)
    sim.run(t[-1], obs)

    return Result(w, obs)
