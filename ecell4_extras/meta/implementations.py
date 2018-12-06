# coding: utf-8

import math
import collections
from coordinator import SimulatorAdapter, DiscreteEvent, DiscreteTimeEvent, EventKind
from ecell4.core import Real3, Integer3, Species, ReactionRule, GSLRandomNumberGenerator, ParticleID, Voxel, Particle
from ecell4 import gillespie, meso, spatiocyte, egfrd, ode


def simulator_event(sim):
    if isinstance(sim, ode.ODESimulator):
        return ODEEvent(sim)
    elif isinstance(sim, gillespie.GillespieSimulator):
        return GillespieEvent(sim)
    elif isinstance(sim, meso.MesoscopicSimulator):
        return MesoscopicEvent(sim)
    elif isinstance(sim, spatiocyte.SpatiocyteSimulator):
        return SpatiocyteEvent(sim)
    elif isinstance(sim, egfrd.EGFRDSimulator):
        return EGFRDEvent(sim)
    raise ValueError("{} not supported.".format(repr(sim)))

class ODESimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs, last_reactions=[], rng=GSLRandomNumberGenerator()):
        SimulatorAdapter.__init__(self, lhs, rhs)
        self.rng = rng
        self.__last_reactions = last_reactions
        assert isinstance(self.lhs.sim, ode.ODESimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return self.__last_reactions
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                return meso.ReactionInfo(
                    ri.t(), ri.reactants(), ri.products(),
                    self.rng.uniform_int(
                        0, self.rhs.world.num_subvolumes() - 1))
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                coord = self.rng.uniform_int(0, self.rhs.world.size() - 1)
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                lengths = self.lhs.world.edge_lengths()
                pos = Real3(self.rng.uniform(0, 1) * lengths[0],
                            self.rng.uniform(0, 1) * lengths[1],
                            self.rng.uniform(0, 1) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        raise ValueError("Not supported yet.")

class ODEEvent(DiscreteTimeEvent):

    def __init__(self, sim, dt=None):
        if dt is None:
            dt = sim.dt()
        DiscreteTimeEvent.__init__(self, sim, dt)
        self.__last_reactions = []

    def step(self):
        DiscreteTimeEvent.step(self)
        self.__last_reactions = self.generate_reactions()
        self.event_kind = EventKind.REACTION_EVENT
        # if len(self.__last_reactions) > 0:
        #     self.event_kind = EventKind.REACTION_EVENT
        # else:
        #     self.event_kind = EventKind.STEP_EVENT

    def sync(self):
        self.__last_reactions = []
        self.event_kind = EventKind.NO_EVENT

        dirty = False
        for sp in self.world.list_species():
            if self.own(sp):
                continue
            value = self.world.get_value_exact(sp)
            if value >= 1.0:
                self.world.set_value(sp, value - math.floor(value))
                dirty = True

        if dirty:
            self.sim.initialize()

    def apply(self, t, ri):
        reactants = [sp for sp in ri.reactants() if self.own(sp)]
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(reactants) == 0 and len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in reactants:
            self.world.remove_molecules(sp, 1)
        for sp in products:
            self.world.add_molecules(sp, 1)
        return True

    def adapter(self, rhs):
        assert self.sim != rhs.sim
        return ODESimulatorAdapter(self, rhs, self.__last_reactions)

    def generate_reactions(self):
        retval = []
        for sp in self.world.list_species():
            if self.own(sp):
                continue
            value = self.world.get_value_exact(sp)
            if value >= 1.0:
                rr = ReactionRule([], [sp], 0.0)
                ri = gillespie.ReactionInfo(self.sim.t(), [], [sp])
                for i in range(int(value)):
                    retval.append((rr, ri))
        return retval

class GillespieSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, gillespie.GillespieSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = ri.reactants()
                if len(reactants) < 2:
                    return ri
                else:
                    assert len(reactants) == 2
                    sp = self.lhs.owe(reactants[1])
                    if sp is None or not self.rhs.own(sp):
                        return ri
                    else:
                        reactants[1] = sp
                        return gillespie.ReactionInfo(ri.t(), reactants, ri.products())
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = ri.reactants()
                if len(reactants) < 2:
                    coord = self.lhs.world.rng().uniform_int(
                        0, self.rhs.world.num_subvolumes() - 1)
                else:
                    assert len(reactants) == 2
                    sp = self.lhs.owe(reactants[1])
                    if sp is None or not self.rhs.own(sp):
                        coord = self.lhs.world.rng().uniform_int(
                            0, self.rhs.world.num_subvolumes() - 1)
                    else:
                        reactants[1] = sp
                        coords = self.rhs.world.list_coordinates_exact(sp)
                        coord = coords[self.lhs.world.rng().uniform_int(0, len(coords) - 1)]
                return meso.ReactionInfo(ri.t(), reactants, ri.products(), coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                if len(ri.reactants()) < 2:
                    coord = self.lhs.world.rng().uniform_int(0, self.rhs.world.size() - 1)
                    reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                 for sp in ri.reactants()]
                    products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                 for sp in ri.products()]
                else:
                    assert len(ri.reactants()) == 2
                    sp = self.lhs.owe(ri.reactants()[1])
                    if sp is None or not self.rhs.own(sp):
                        coord = self.lhs.world.rng().uniform_int(0, self.rhs.world.size() - 1)
                        reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                     for sp in ri.reactants()]
                        products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                     for sp in ri.products()]
                    else:
                        voxels = self.rhs.world.list_voxels_exact(sp)
                        v = voxels[self.lhs.world.rng().uniform_int(0, len(voxels) - 1)]
                        reactants = [(ParticleID(), Voxel(ri.reactants()[0], v[1].coordinate(), 0, 0)), v]
                        products = [(ParticleID(), Voxel(sp, v[1].coordinate(), 0, 0))
                                     for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                rng = self.lhs.world.rng()
                lengths = self.lhs.world.edge_lengths()
                if len(ri.reactants()) < 2:
                    pos = Real3(rng.uniform(0, 1) * lengths[0],
                                rng.uniform(0, 1) * lengths[1],
                                rng.uniform(0, 1) * lengths[2])
                    reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                                 for sp in ri.reactants()]
                    products = [(ParticleID(), Particle(sp, pos, 0, 0))
                                 for sp in ri.products()]
                else:
                    assert len(ri.reactants()) == 2
                    sp = self.lhs.owe(ri.reactants()[1])
                    if sp is None or not self.rhs.own(sp):
                        pos = Real3(rng.uniform(0, 1) * lengths[0],
                                    rng.uniform(0, 1) * lengths[1],
                                    rng.uniform(0, 1) * lengths[2])
                        reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                                     for sp in ri.reactants()]
                        products = [(ParticleID(), Particle(sp, pos, 0, 0))
                                     for sp in ri.products()]
                    else:
                        particles = self.rhs.world.list_particles_exact(sp)
                        p = particles[rng.uniform_int(0, len(particles) - 1)]
                        reactants = [(ParticleID(), Particle(sp, p[1].position(), 0, 0)), p]
                        products = [(ParticleID(), Particle(sp, p[1].position(), 0, 0))
                                     for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet.")

class GillespieEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for sp in ri.products():
            if self.own(sp):
                continue
            dirty = True
            self.world.remove_molecules(sp, 1)

        if dirty:
            self.sim.initialize()

    def apply(self, t, ri):
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.world.add_molecules(sp, 1)
        return True

    def mirror(self, t, interrupter, src, dst):
        value1 = interrupter.world.get_value_exact(src)
        assert value1 >= 0
        value2 = self.world.get_value_exact(dst)
        if value1 != value2:
            self.sim.step(t)
            self.world.set_value(dst, value1)
            return True
        return False

    def adapter(self, rhs):
        assert self.sim != rhs.sim
        return GillespieSimulatorAdapter(self, rhs)

class MesoscopicWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def list_coordinates_exact(self, sp):
        assert isinstance(self.rhs.world, meso.MesoscopicWorld)
        ratios = [
            float(y) / float(x)
            for x, y in zip(self.lhs.world.matrix_sizes(), self.rhs.world.matrix_sizes())]
        rng = self.lhs.world.rng()
        def convert(coord):
            g1 = self.lhs.world.coord2global(coord)
            g2 = Integer3(*[int(rng.uniform(g1[i] * ratios[i], (g1[i] + 1) * ratios[i])) for i in range(3)])
            return self.rhs.world.global2coord(g2)
        coords = [convert(coord) for coord in self.lhs.world.list_coordinates_exact(sp)]
        coords.sort()
        return coords

    def __getattr__(self, name):
        return getattr(self.lhs.world, name)

class MesoscopicSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, meso.MesoscopicSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            ratios = [
                float(y) / float(x)
                for x, y in zip(self.lhs.world.matrix_sizes(), self.rhs.world.matrix_sizes())]
            rng = self.lhs.world.rng()
            def convert(ri):
                reactants = ri.reactants()
                coord = ri.coordinate()
                g1 = self.lhs.world.coord2global(coord)
                if len(reactants) < 2:
                    g2 = Integer3(*[int(rng.uniform(g1[i] * ratios[i], (g1[i] + 1) * ratios[i])) for i in range(3)])
                    newcoord = self.rhs.world.global2coord(g2)
                else:
                    assert len(reactants) == 2
                    sp = self.lhs.owe(reactants[1])
                    if sp is None or not self.rhs.own(sp):
                        g2 = Integer3(*[int(rng.uniform(g1[i] * ratios[i], (g1[i] + 1) * ratios[i])) for i in range(3)])
                        newcoord = self.rhs.world.global2coord(g2)
                    else:
                        coords = [c for c in self.rhs.world.list_coordinates_exact(sp) if all(g1[i] * ratios[i] <= self.rhs.world.coord2global(c)[i] < (g1[i] + 1) * ratios[i] for i in range(3))]
                        assert len(coords) > 0
                        newcoord = coords[rng.uniform_int(0, len(coords) - 1)]
                        reactants[1] = sp
                        # print(tuple(self.lhs.world.coord2global(coord)), tuple(self.rhs.world.coord2global(newcoord)))
                return meso.ReactionInfo(ri.t(), reactants, ri.products(), newcoord)
            return [(rr, convert(ri))
                    for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return [(rr, gillespie.ReactionInfo(ri.t(), ri.reactants(), ri.products()))
                    for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                reactants = ri.reactants()
                g = self.lhs.world.coord2global(ri.coordinate())
                rng = self.lhs.world.rng()
                lengths = self.lhs.world.subvolume_edge_lengths()
                if len(reactants) < 2:
                    pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                                (g[1] + rng.uniform(0, 1)) * lengths[1],
                                (g[2] + rng.uniform(0, 1)) * lengths[2])
                    assert ri.coordinate() == self.lhs.world.position2coordinate(pos)
                    coord = self.rhs.world.position2coordinate(pos)  #XXX: This may cause an overlap
                    # assert ri.coordinate() == self.lhs.world.position2coordinate(self.rhs.world.coordinate2position(coord))  #XXX: This is not always True.
                    reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                 for sp in ri.reactants()]
                    products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                 for sp in ri.products()]
                else:
                    assert len(reactants) == 2
                    sp = self.lhs.owe(reactants[1])
                    if sp is None or not self.rhs.own(sp):
                        pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                                    (g[1] + rng.uniform(0, 1)) * lengths[1],
                                    (g[2] + rng.uniform(0, 1)) * lengths[2])
                        assert ri.coordinate() == self.lhs.world.position2coordinate(pos)
                        coord = self.rhs.world.position2coordinate(pos)  #XXX: This may cause an overlap
                        # assert ri.coordinate() == self.lhs.world.position2coordinate(self.rhs.world.coordinate2position(coord))  #XXX: This is not always True.
                        reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                     for sp in ri.reactants()]
                        products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                                     for sp in ri.products()]
                    else:
                        voxels = [(pid, v) for pid, v in self.rhs.world.list_voxels_exact(sp) if self.lhs.world.position2coordinate(self.rhs.world.coordinate2position(v.coordinate())) == ri.coordinate()]
                        v = voxels[rng.uniform_int(0, len(voxels) - 1)]
                        reactants = [(ParticleID(), Voxel(ri.reactants()[0], v[1].coordinate(), 0, 0)), v]
                        products = [(ParticleID(), Voxel(sp, v[1].coordinate(), 0, 0))
                                     for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                reactants = ri.reactants()
                g = self.lhs.world.coord2global(ri.coordinate())
                rng = self.lhs.world.rng()
                lengths = self.lhs.world.subvolume_edge_lengths()
                if len(reactants) < 2:
                    pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                                (g[1] + rng.uniform(0, 1)) * lengths[1],
                                (g[2] + rng.uniform(0, 1)) * lengths[2])
                    reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                                 for sp in ri.reactants()]
                    products = [(ParticleID(), Particle(sp, pos, 0, 0))
                                 for sp in ri.products()]
                else:
                    assert len(reactants) == 2
                    sp = self.lhs.owe(reactants[1])
                    if sp is None or not self.rhs.own(sp):
                        pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                                    (g[1] + rng.uniform(0, 1)) * lengths[1],
                                    (g[2] + rng.uniform(0, 1)) * lengths[2])
                        reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                                     for sp in ri.reactants()]
                        products = [(ParticleID(), Particle(sp, pos, 0, 0))
                                     for sp in ri.products()]
                    else:
                        particles = [(pid, p) for pid, p in self.rhs.world.list_particles_exact(sp) if self.lhs.world.position2coordinate(p.position()) == ri.coordinate()]
                        p = particles[rng.uniform_int(0, len(particles) - 1)]
                        reactants = [(ParticleID(), Particle(sp, p[1].position(), 0, 0)), p]
                        products = [(ParticleID(), Particle(sp, p[1].position(), 0, 0))
                                     for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

    def world(self):
        return MesoscopicWorldAdapter(self.lhs, self.rhs)

class MesoscopicEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]
        coord = ri.coordinate()

        dirty = False
        for sp in ri.products():
            if self.own(sp):
                continue
            dirty = True
            self.world.remove_molecules(sp, 1, coord)

        if dirty:
            self.sim.initialize()

    def apply(self, t, ri):
        reactants = [sp for sp in ri.reactants() if self.own(sp)]
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(reactants) == 0 and len(products) == 0:
            return False
        coord = ri.coordinate()
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in reactants:
            self.world.remove_molecules(sp, 1, coord)
        for sp in products:
            self.world.add_molecules(sp, 1, coord)
        return True

    def mirror(self, t, interrupter, src, dst):
        coords1 = self.world.list_coordinates_exact(dst)
        coords2 = interrupter(self).world().list_coordinates_exact(src)
        if coords1 != coords2:
            self.sim.step(t)
            for c1, c2 in zip(coords1, coords2):
                if c1 != c2:
                    self.world.remove_molecules(dst, 1, c1)
                    self.world.add_molecules(dst, 1, c2)
            if len(coords1) > len(coords2):
                for c1 in coords1[len(coords2):]:
                    self.world.remove_molecules(dst, 1, c1)
            elif len(coords1) < len(coords2):
                for c2 in coords2[len(coords1):]:
                    self.world.add_molecules(dst, 1, c2)
            return True
        return False

    def adapter(self, rhs):
        assert self.sim != rhs.sim
        return MesoscopicSimulatorAdapter(self, rhs)

class SpatiocyteWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def list_coordinates_exact(self, sp):
        assert isinstance(self.rhs.world, meso.MesoscopicWorld)
        coords = [self.rhs.world.position2coordinate(
                      self.lhs.world.coordinate2position(v.coordinate()))
                  for pid, v in self.lhs.world.list_voxels_exact(sp)]
        coords.sort()
        return coords

    def __getattr__(self, name):
        return getattr(self.lhs.world, name)

class SpatiocyteSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, spatiocyte.SpatiocyteSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            raise RuntimeError("Not supported yet.")
            # return self.lhs.sim.last_reactions()
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                assert len(products) > 0
                pos = self.lhs.world.coordinate2position(ri.products()[0][1].coordinate())
                coord = self.rhs.world.position2coordinate(pos)
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                reactants = [(pid, Particle(v.species(), self.lhs.world.coordinate2position(v.coordinate()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Particle(v.species(), self.lhs.world.coordinate2position(v.coordinate()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

    def world(self):
        return SpatiocyteWorldAdapter(self.lhs, self.rhs)

class SpatiocyteEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for pid, v in ri.products():
            if self.own(v.species()):
                continue
            dirty = True
            self.world.remove_voxel(pid)

        if dirty:
            self.sim.initialize()

    def apply(self, t, ri):
        reactants = [(pid, v) for pid, v in ri.reactants() if self.own(v.species())]
        products = [(pid, v) for pid, v in ri.products() if self.own(v.species())]
        if len(reactants) == 0 and len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, v in reactants:
            self.world.remove_voxel(pid)
        for pid, v in products:
            self.world.new_voxel(v.species(), v.coordinate())
        return True

    def adapter(self, rhs):
        assert self.sim != rhs.sim
        return SpatiocyteSimulatorAdapter(self, rhs)

class EGFRDWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def list_coordinates_exact(self, sp):
        assert isinstance(self.rhs.world, meso.MesoscopicWorld)
        coords = [self.rhs.world.position2coordinate(p.position())
                  for pid, p in self.lhs.world.list_particles_exact(sp)]
        coords.sort()
        return coords

    def __getattr__(self, name):
        return getattr(self.lhs.world, name)

class EGFRDSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, egfrd.EGFRDSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                reactants = [(pid, Voxel(v.species(), self.rhs.world.position2coordinate(v.position()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Voxel(v.species(), self.rhs.world.position2coordinate(v.position()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                assert len(products) > 0
                coord = self.rhs.world.position2coordinate(ri.products()[0][1].position())
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            return self.lhs.sim.last_reactions()
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

    def world(self):
        return EGFRDWorldAdapter(self.lhs, self.rhs)

class EGFRDEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for pid, p in ri.products():
            if self.own(p.species()):
                continue
            dirty = True
            self.world.remove_particle(pid)

        if dirty:
            self.sim.initialize()

    def apply(self, t, ri):
        reactants = [(pid, p) for pid, p in ri.reactants() if self.own(p.species())]
        products = [(pid, p) for pid, p in ri.products() if self.own(p.species())]
        if len(reactants) == 0 and len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, p in reactants:
            self.world.remove_particle(pid)
        for pid, p in products:
            self.world.new_particle(p.species(), p.position())
        return True

    def adapter(self, rhs):
        assert self.sim != rhs.sim
        return EGFRDSimulatorAdapter(self, rhs)
