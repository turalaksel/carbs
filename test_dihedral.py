from hoomd import *
from hoomd import md
import numpy as np

# Start HOOMD
context.initialize("");

num_part = 10

snapshot = data.make_snapshot(N = num_part,
                              box = data.boxdim(Lx=40, Ly=40, Lz=40),
                              particle_types=['C','A'],
                              bond_types = ['polymer'],
                              angle_types = [],
                              dihedral_types = ['dihedral'],
                              improper_types = []);

#positions
part_pos = [[0,0,i] for i in range(0,num_part)]
snapshot.particles.position[:] = part_pos
#intertia
snapshot.particles.moment_inertia[:] = [[10,10,10]]*num_part
#type
snapshot.particles.typeid[:] = [0];
#bonds
bonds = [[n, min(n+1, num_part)] for n in range(0, num_part - 1, 1)]
snapshot.bonds.resize(num_part - 1);
snapshot.bonds.group[:] = bonds

system = init.read_snapshot(snapshot);

nl = md.nlist.cell();

rigid = md.constrain.rigid();

rigid.set_param('C',
                types=['A']*2,
                positions = [[0.5,0,0],[1.0,0,0]]
                           );

rigid.create_bodies()

# nl.reset_exclusions(exclusions = []);

harmonic = md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=50.0, r0=1);

snapshot.dihedrals.resize(num_part);
for i in range(num_part):
    snapshot.dihedrals.group[i,:] = [num_part + 2*i, i, i+1, num_part + 2*(i+1)];

system.restore_snapshot(snapshot)

# LJ interactions
lj = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
lj.set_params(mode='shift')
lj.pair_coeff.set(['C', 'A'], ['C', 'A'], epsilon=1.0, sigma=1.0)

md.integrate.mode_standard(dt=0.005, aniso=True);

rigid = group.rigid_center();
md.integrate.langevin(group=rigid, kT=0.1, seed=42);

hoomd.dump.gsd("dihedral.gsd",
               period=1e3,
               group=hoomd.group.all(),
               overwrite=True);
hoomd.run(10e4);
