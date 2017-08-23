from hoomd import *
from hoomd import md
import numpy as np
import random
import particlesFromPDB as fromPDB

# Start HOOMD
context.initialize("");

num_part = 42 #even !!

snapshot = data.make_snapshot(N = num_part,
                              box = data.boxdim(Lx=60, Ly=60, Lz=40),
                              particle_types=['C','A'],
                              bond_types = ['polymer'],
                              angle_types = ['angle'],
                              dihedral_types = ['dihedral'],
                              improper_types = []);

#positions
def ran(i):
    return(random.uniform(-0.2, 0.2))


part_pos = [[ran(i),ran(i),i] for i in range(-int(num_part/2),int(num_part/2))]
snapshot.particles.position[:] = part_pos
#inertia
snapshot.particles.moment_inertia[:] = [[10,10,10]]*num_part
#type
snapshot.particles.typeid[:] = [0];
#bonds
bonds = [[n, min(n+1, num_part)] for n in range(0, num_part - 1, 1)]

snapshot.bonds.resize(num_part - 1);
snapshot.bonds.group[:] = bonds

#read the snapshot and create neighbor list
system = init.read_snapshot(snapshot);
nl = md.nlist.cell();

#rigid
rigid = md.constrain.rigid();
rigid.set_param('C', types=['A'], positions = [[1,0,0]]);
rigid.create_bodies()

#fene / harmonic
harmonic1 = md.bond.harmonic()
harmonic1.bond_coeff.set('polymer', k=25.0, r0=1.0)

#angle for backbone
for i in range(num_part-2):
    system.angles.add('angle', i, i+1, i+2)
harmonic2 = md.angle.harmonic()
harmonic2.angle_coeff.set('angle', k=30.0, t0=2.9)

#dihedral
def harmonicAngle(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)
for i in range(num_part-4):
    system.dihedrals.add('dihedral', i,i+1,i+2,i+3)
dtable = md.dihedral.table(width=1000)
dtable.dihedral_coeff.set('dihedral', func=harmonicAngle, coeff=dict(kappa=15, theta0=0.2))


# system.particles[0].position

# LJ interactions
wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
wca.set_params(mode='shift')
wca.pair_coeff.set(['C', 'A'],['C', 'A'], epsilon=1.0, sigma=1.0, r_cut=0.7**(1/6))

md.integrate.mode_standard(dt=0.003, aniso=True);

rigid = group.rigid_center();
md.integrate.langevin(group=rigid, kT=0.1, seed=42);

dump.gsd("dihedral.gsd",
               period=1e3,
               group=group.all(),
               overwrite=True);
run(10e4);
