# -*- coding: iso-8859-1 -*-
# Maintainer: joaander

from hoomd import *
from hoomd import md
import numpy

context.initialize("");

snap = data.make_snapshot(N=40,
                          box=data.boxdim(L=100),
                          particle_types = ['A'],
                          bond_types = [],
                          angle_types = [],
                          dihedral_types = ['DihA'],
                          improper_types = [])

snap.dihedrals.resize(10);

for i in range(10):
    x = numpy.array([i, 0, 0], dtype=numpy.float32)
    snap.particles.position[4*i+0,:] = x;
    x += numpy.random.random(3)
    snap.particles.position[4*i+1,:] = x;
    x += numpy.random.random(3)
    snap.particles.position[4*i+2,:] = x;
    x += numpy.random.random(3)
    snap.particles.position[4*i+3,:] = x;

    snap.dihedrals.group[i,:] = [4*i+0, 4*i+1, 4*i+2, 4*i+3];

init.read_snapshot(snap)

# context.current.sorter.set_params(grid=8)

harmonic = md.dihedral.harmonic();
harmonic.dihedral_coeff.set('DihA', k=1.0, d=1, n=4)

all = group.all();
md.integrate.mode_standard(dt=0.005);
md.integrate.nve(all);

dump.gsd("dihedral.gsd",
               period=1e3,
               group=all,
               overwrite=True);

run(1e6);
