from __future__ import print_function

import ecell4_base

from .viz import plot_number_observer, plot_trajectory, plot_world
from .simulation import load_world


def dump_model(m):
    res = [sp.serial() + "|" + str(dict(sp.list_attributes())) for sp in m.species_attributes()]
    res += [rr.as_string() for rr in m.reaction_rules()]
    # return res
    print('\n'.join(res))

def show(target, *args, **kwargs):
    """
    An utility function to display the given target object in the proper way.

    Paramters
    ---------
    target : NumberObserver, TrajectoryObserver, World, str
        When a NumberObserver object is given, show it with viz.plot_number_observer.
        When a TrajectoryObserver object is given, show it with viz.plot_trajectory_observer.
        When a World or a filename suggesting HDF5 is given, show it with viz.plot_world.

    """
    if isinstance(target, (ecell4_base.core.FixedIntervalNumberObserver, ecell4_base.core.NumberObserver, ecell4_base.core.TimingNumberObserver, )):
        plot_number_observer(target, *args, **kwargs)
    elif isinstance(target, (ecell4_base.core.FixedIntervalTrajectoryObserver, ecell4_base.core.FixedIntervalTrackingObserver)):
        plot_trajectory(target, *args, **kwargs)
    elif isinstance(target, (ecell4_base.ode.ODEWorld, ecell4_base.gillespie.GillespieWorld, ecell4_base.spatiocyte.SpatiocyteWorld, ecell4_base.meso.MesoscopicWorld, ecell4_base.bd.BDWorld, ecell4_base.egfrd.EGFRDWorld)):
        plot_world(target, *args, **kwargs)
    elif isinstance(target, (ecell4_base.core.Model, ecell4_base.core.NetworkModel, ecell4_base.core.NetfreeModel)):
        dump_model(target)
    elif isinstance(target, str):
        try:
            w = simulation.load_world(target)
        except RuntimeError as e:
            raise ValueError("The given target [{}] is not supported.".format(repr(target)))
        else:
            show(w, *args, **kwargs)
    else:
        raise ValueError("The given target [{}] is not supported.".format(repr(target)))
