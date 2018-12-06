# coding: utf-8

import numpy

from ecell4.core import Real3, Integer3, Species
from ecell4 import gillespie, meso, spatiocyte, egfrd, ode
from ecell4.util import species_attributes, reaction_rules, get_model

from coordinator import Coordinator, EventKind
from implementations import simulator_event, ODEEvent
from logger import Logger


def test1():
    D, radius = 1, 0.005
    edge_lengths = Real3(1, 1, 1)

    with species_attributes():
        A1 | A2 | B1 | B2 | C1 | C2 | D1 | D2 | E1 | E2 | {
            "D": str(D), "radius": str(radius)}

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        B1 == B2 | (1.0, 1.0)
        C1 == C2 | (1.0, 1.0)
        D1 == D2 | (1.0, 1.0)
        E1 == E2 | (1.0, 1.0)

        A1 == B1 | (1.0, 1.0)
        A1 == C1 | (1.0, 1.0)
        A1 == D1 | (1.0, 1.0)
        A1 == E1 | (1.0, 1.0)
        B1 == C1 | (1.0, 1.0)
        B1 == D1 | (1.0, 1.0)
        B1 == E1 | (1.0, 1.0)
        C1 == D1 | (1.0, 1.0)
        C1 == E1 | (1.0, 1.0)
        D1 == E1 | (1.0, 1.0)


    m = get_model()

    w1 = gillespie.GillespieWorld(edge_lengths)
    w1.bind_to(m)
    sim1 = gillespie.GillespieSimulator(w1)

    w2 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
    w2.bind_to(m)
    sim2 = meso.MesoscopicSimulator(w2)

    w3 = spatiocyte.SpatiocyteWorld(edge_lengths, radius)
    w3.bind_to(m)
    sim3 = spatiocyte.SpatiocyteSimulator(w3)

    w4 = egfrd.EGFRDWorld(edge_lengths, Integer3(4, 4, 4))
    w4.bind_to(m)
    sim4 = egfrd.EGFRDSimulator(m, w4)

    w5 = ode.ODEWorld(edge_lengths)
    w5.bind_to(m)
    sim5 = ode.ODESimulator(w5)
    sim5.set_dt(0.01)

    owner = Coordinator()
    owner.add_event(simulator_event(sim1)).add(('A1', 'A2'))
    owner.add_event(simulator_event(sim2)).add(('B1', 'B2'))
    owner.add_event(simulator_event(sim3)).add(('C1', 'C2'))
    owner.add_event(simulator_event(sim4)).add(('D1', 'D2'))
    owner.add_event(simulator_event(sim5)).add(('E1', 'E2'))

    owner.set_value(Species("A1"), 300)
    # owner.set_value(Species("A1"), 60)
    # owner.set_value(Species("B1"), 60)
    # owner.set_value(Species("C1"), 60)
    # owner.set_value(Species("D1"), 60)
    # owner.set_value(Species("E1"), 60)

    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "E1", "E2"))

    logger.log()
    while owner.step(3):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()
    logger.log()

    logger.savefig()
    logger.savetxt()

def test2():
    edge_lengths = Real3(1, 1, 1)

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        E1 == E2 | (1.0, 1.0)
        A1 == E1 | (1.0, 1.0)

    m = get_model()

    w1 = gillespie.GillespieWorld(edge_lengths)
    w1.bind_to(m)
    sim1 = gillespie.GillespieSimulator(w1)

    w2 = ode.ODEWorld(edge_lengths)
    w2.bind_to(m)
    sim2 = ode.ODESimulator(w2)

    owner = Coordinator()
    owner.add_event(simulator_event(sim1)).add(('A1', 'A2'))
    owner.add_event(ODEEvent(sim2, 0.01)).add(('E1', 'E2'))
    owner.set_value(Species("A1"), 120)
    # owner.set_value(Species("A1"), 60)
    # owner.set_value(Species("E1"), 60)
    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "E1", "E2"))
    logger.add("A1", w2)

    logger.log()
    while owner.step(50):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()

    logger.savefig()
    logger.savetxt()

def test3():
    D, radius = 1, 0.005
    edge_lengths = Real3(1, 1, 1)

    with species_attributes():
        A1 | A2 | B1 | B2 | B3 | {
            "D": str(D), "radius": str(radius)}

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        B1 == B2 | (1.0, 1.0)

        A2 + B2_ > B3 | 1.0 / 30.0
        B3 > A2 + B2 | 1.0

    m = get_model()

    w1 = gillespie.GillespieWorld(edge_lengths)
    w1.bind_to(m)
    sim1 = gillespie.GillespieSimulator(w1)

    # w2 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
    # w2.bind_to(m)
    # sim2 = meso.MesoscopicSimulator(w2)
    # w2 = spatiocyte.SpatiocyteWorld(edge_lengths, radius)
    # w2.bind_to(m)
    # sim2 = spatiocyte.SpatiocyteSimulator(w2)
    # w2 = egfrd.EGFRDWorld(edge_lengths, Integer3(4, 4, 4))
    # w2.bind_to(m)
    # sim2 = egfrd.EGFRDSimulator(w2)
    w2 = ode.ODEWorld(edge_lengths)
    w2.bind_to(m)
    sim2 = ode.ODESimulator(w2)
    sim2.set_dt(0.01)

    owner = Coordinator()
    ev1 = simulator_event(sim1)
    ev1.add(('A1', 'A2'))
    ev1.borrow('B2', 'B2_')
    owner.add_event(ev1)
    owner.add_event(simulator_event(sim2)).add(('B1', 'B2', 'B3'))
    owner.set_value(Species("A1"), 60)
    owner.set_value(Species("B1"), 60)
    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "B1", "B2", "B3"))
    logger.add("B2_", w1)

    logger.log()
    while owner.step(10):
    # while owner.step(50):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()
    logger.log()

    logger.savefig()
    logger.savetxt()

def test4():
    D, radius = 1, 0.005
    edge_lengths = Real3(1, 1, 1)

    with species_attributes():
        A1 | A2 | B1 | B2 | {"D": str(D), "radius": str(radius)}

    with reaction_rules():
        A1 + B1_ > B2 | 0.04483455086786913 > B1 + A2 | 1.35
        B2 > B1 + A1 | 1.5

    m = get_model()

    w1 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
    w1.bind_to(m)
    sim1 = meso.MesoscopicSimulator(w1)

    # w2 = spatiocyte.SpatiocyteWorld(edge_lengths, radius)
    # w2.bind_to(m)
    # sim2 = spatiocyte.SpatiocyteSimulator(w2)
    w2 = egfrd.EGFRDWorld(edge_lengths, Integer3(4, 4, 4))
    w2.bind_to(m)
    sim2 = egfrd.EGFRDSimulator(w2)

    w1.add_molecules(Species("A1"), 120)
    w2.add_molecules(Species("B1"), 60)

    owner = Coordinator()
    ev1 = simulator_event(sim1)
    ev1.add(('A1', 'A2'))
    ev1.borrow('B1', 'B1_')
    owner.add_event(ev1)
    owner.add_event(simulator_event(sim2)).add(('B1', 'B2'))
    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "B1", "B2"))
    logger.add("B1_", w1)

    logger.log()
    while owner.step(1):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()
    logger.log()

    logger.savefig()
    logger.savetxt()

def test4():
    D, radius = 1, 0.005
    edge_lengths = Real3(1, 1, 1)

    with species_attributes():
        A1 | A2 | B1 | B2 | B3 | C1 | C2 | C3 | {"D": str(D), "radius": str(radius)}

    with reaction_rules():
        A1 + B1_ > B2 | 0.04483455086786913 > A1 + B1 | 1.35
        B2 > A2 + B1 | 1.5
        A2 + B1_ > B3 | 0.09299017957780264 > A2 + B1 | 1.73
        B3 > A3 + B1 | 15.0

        A3 + C1 > C2 | 0.04483455086786913 > A3 + C1 | 1.35
        C2 > A2 + C1 | 1.5
        A2 + C1 > C3 | 0.09299017957780264 > A2 + C1 | 1.73
        C3 > A1 + C1 | 15.0

    m = get_model()

    w1 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
    w1.bind_to(m)
    sim1 = meso.MesoscopicSimulator(w1)

    w2 = spatiocyte.SpatiocyteWorld(edge_lengths, radius)
    w2.bind_to(m)
    sim2 = spatiocyte.SpatiocyteSimulator(w2)
    # w2 = egfrd.EGFRDWorld(edge_lengths, Integer3(4, 4, 4))
    # w2.bind_to(m)
    # sim2 = egfrd.EGFRDSimulator(w2)

    owner = Coordinator()
    ev1 = simulator_event(sim1)
    ev1.add(('A1', 'A2', 'A3'))
    ev1.add(('C1', 'C2', 'C3'))
    ev1.borrow('B1', 'B1_')
    owner.add_event(ev1)
    owner.add_event(simulator_event(sim2)).add(('B1', 'B2', 'B3'))
    owner.set_value(Species("A1"), 120)
    owner.set_value(Species("B1"), 30)
    owner.set_value(Species("C1"), 30)
    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"))
    logger.add("B1_", w1)

    logger.log()
    while owner.step(1):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()
    logger.log()

    logger.savefig()
    logger.savetxt()

def test5():
    edge_lengths = Real3(1, 1, 1)

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        B1 == B2 | (1.0, 1.0)
        A1 == B1 | (1.0, 1.0)

        A2 + B2_ > C1 | 1.0 / 30.0
        C1 > A2 + B2 | 1.0

    m = get_model()

    w1 = meso.MesoscopicWorld(edge_lengths, Integer3(3, 3, 3))
    w1.bind_to(m)
    sim1 = meso.MesoscopicSimulator(w1)

    w2 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
    w2.bind_to(m)
    sim2 = meso.MesoscopicSimulator(w2)

    owner = Coordinator()
    ev1 = simulator_event(sim1)
    ev1.add(('A1', 'A2', 'C1'))
    ev1.borrow('B2', 'B2_')
    owner.add_event(ev1)
    owner.add_event(simulator_event(sim2)).add(('B1', 'B2'))
    owner.set_value(Species("A1"), 120)
    # owner.set_value(Species("A1"), 60)
    # owner.set_value(Species("B1"), 60)
    owner.initialize()

    logger = Logger(owner, ("A1", "A2", "C1", "B1", "B2"))
    logger.add('B2_', w1)

    logger.log()
    while owner.step(5):
        if owner.last_event.event_kind == EventKind.REACTION_EVENT:
            logger.log()

    logger.savefig()
    logger.savetxt()


if __name__ == "__main__":
    # test1()
    # test2()
    # test3()
    # test4()
    test5()
