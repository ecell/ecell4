from __future__ import print_function
import sys

import ecell4_base

from .viz import plot_number_observer, plot_trajectory, plot_world
from .simulation import load_world


def _as_string(val):
    if isinstance(val, ecell4_base.core.Quantity_Integer):
        if val.units == '':
            return str(val.magnitude)
        return "Quantity_Integer({:d}, '{:s}')".format(val.magnitude, val.units)
    elif isinstance(val, ecell4_base.core.Quantity_Real):
        if val.units == '':
            return str(val.magnitude)
        return "Quantity_Real({:g}, '{:s}')".format(val.magnitude, val.units)
    elif isinstance(val, str):
        return repr(val)
    return str(val)

def _show_species(sp, f=sys.stdout):
    attributes = ', '.join("'{:s}': {:s}".format(key, _as_string(val)) for key, val in sp.list_attributes())
    f.write("{:s} | {{{:s}}}\n".format(sp.serial(), attributes))

def _show_reaction_rule(rr, f=sys.stdout):
    if not rr.has_descriptor():
        k = rr.get_k()
        lhs = ' + '.join(sp.serial() for sp in rr.reactants())
        rhs = ' + '.join(sp.serial() for sp in rr.products())
        k = _as_string(k)
    else:
        desc = rr.get_descriptor()
        k = rr.get_k()
        lhs = ' + '.join('{:g} * {:s}'.format(coef, sp.serial()) if coef != 1 else sp.serial() for sp, coef in zip(rr.reactants(), desc.reactant_coefficients()))
        rhs = ' + '.join('{:g} * {:s}'.format(coef, sp.serial()) if coef != 1 else sp.serial() for sp, coef in zip(rr.products(), desc.product_coefficients()))
        k = desc.as_string()

    f.write("{:s} > {:s} | {:s}\n".format(lhs, rhs, k))

def _show_model(m, f=sys.stdout):
    for sp in m.species_attributes():
        _show_species(sp, f)
    f.write('\n')
    for rr in m.reaction_rules():
        _show_reaction_rule(rr, f)

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
        _show_model(target)
    elif isinstance(target, str):
        try:
            w = simulation.load_world(target)
        except RuntimeError as e:
            raise ValueError("The given target [{}] is not supported.".format(repr(target)))
        else:
            show(w, *args, **kwargs)
    else:
        raise ValueError("The given target [{}] is not supported.".format(repr(target)))
