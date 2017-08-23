from hoomd import *
from hoomd import md
import numpy as np
import random

# Start HOOMD
context.initialize("");

num_part = 60 #even !!

snapshot = data.make_snapshot(N = num_part,
                              box = data.boxdim(Lx=num_part, Ly=num_part, Lz=num_part + 30),
                              particle_types=['C','a','b'], #center, sidea, sideb
                              bond_types = ['polymer'],
                              angle_types = ['angle1', 'angle2'],
                              dihedral_types = ['dihedral1','dihedral21','dihedral22','dihedral31','dihedral32'],
                              improper_types = []);
############ INIT ############
#positions
def ran():
    return(round(random.uniform(-0.1, 0.1),3))
part_pos = [[0.,0.,i] for i in range(-int(num_part/2),int(num_part/2))]
snapshot.particles.position[:] = part_pos
#inertia
snapshot.particles.moment_inertia[:] = [[1./3.,1./3.,1./3.]]*num_part
#type
snapshot.particles.typeid[:] = [0];
#bonds
bonds = [[n, min(n+1, num_part)] for n in range(0, num_part - 1, 1)]
snapshot.bonds.resize(num_part - 1);
snapshot.bonds.group[:] = bonds
#read the snapshot and create neighbor list
system = init.read_snapshot(snapshot);
nl = md.nlist.cell();
############ BONDS ############
#rigid
rigid = md.constrain.rigid();
rigid.set_param('C', types=['a','b'], positions = [[-0.1,0.,0.],[0,-1.0,0.]]);
rigid.create_bodies()

# #fene / harmonic
harmonic1 = md.bond.harmonic()
harmonic1.bond_coeff.set('polymer', k=50.0, r0=1.0)

#angle
for i in range(num_part-1):
    a1 = num_part + 2*i
    b1 = num_part + 2*i + 1
    c1 = i
    a2 = num_part + 2*i + 2
    b2 = num_part + 2*i + 3
    c2 = i + 1
    system.angles.add('angle1', c1, c2, a2)
    # system.angles.add('angle1', c2, c1, a1)
    system.angles.add('angle2', c1, c2, b2)
    # system.angles.add('angle2', c2, c1, b1)
harmonic2 = md.angle.harmonic()
harmonic2.angle_coeff.set('angle1', k=00.0, t0=3.1415/2.-0.1) #raise
harmonic2.angle_coeff.set('angle2', k=00.0, t0=3.1415/2.-0.5)

#dihedrals
def harmonicAngle(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)
for i in range(num_part-1):
    a1 = num_part + 2*i
    b1 = num_part + 2*i + 1
    c1 = i
    a2 = num_part + 2*i + 2
    b2 = num_part + 2*i + 3
    c2 = i + 1
    system.dihedrals.add('dihedral1', b1, c1, c2, b2)
    system.dihedrals.add('dihedral21', b1, c1, a1, c2)
    system.dihedrals.add('dihedral22', b2, c2, a2, c1)
    system.dihedrals.add('dihedral31', a1, c1, b1, c2)
    system.dihedrals.add('dihedral32', a2, c2, b2, c1)



dtable = md.dihedral.table(width=1000)
#
# dtable.dihedral_coeff.set('angle1', func=harmonicAngle, coeff=dict(kappa=20, theta0=0.))
# dtable.dihedral_coeff.set('angle2', func=harmonicAngle, coeff=dict(kappa=20, theta0=0.))

dtable.dihedral_coeff.set('dihedral1', func=harmonicAngle, coeff=dict(kappa=30, theta0=-0.34))
dtable.dihedral_coeff.set('dihedral21', func=harmonicAngle, coeff=dict(kappa=30, theta0=1.33))
dtable.dihedral_coeff.set('dihedral22', func=harmonicAngle, coeff=dict(kappa=30, theta0=-1.33))
dtable.dihedral_coeff.set('dihedral31', func=harmonicAngle, coeff=dict(kappa=30, theta0=-1.57))
dtable.dihedral_coeff.set('dihedral32', func=harmonicAngle, coeff=dict(kappa=30, theta0=+1.57))


########## INTERACTIONS ############
# LJ interactions
wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
wca.set_params(mode='shift')
wca.pair_coeff.set('C', 'C', epsilon=1.0, sigma=0.700, r_cut=0.700*2**(1/6))
wca.pair_coeff.set('C', 'b', epsilon=1.0, sigma=0.515, r_cut=0.515*2**(1/6))
wca.pair_coeff.set('b', 'b', epsilon=1.0, sigma=0.330, r_cut=0.330*2**(1/6))
wca.pair_coeff.set('C', 'a', epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
wca.pair_coeff.set('a', 'a', epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
wca.pair_coeff.set('a', 'b', epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
########## INTEGRATION ############

system.particles[0].diameter = 0.1
system.particles[1].diameter = 0.1
system.particles[2].diameter = 0.1
system.particles[3].diameter = 0.1
system.particles[4].diameter = 0.1

md.integrate.mode_standard(dt=0.003, aniso=True);
rigid = group.rigid_center();
md.integrate.langevin(group=rigid, kT=0.3, seed=42);
########## DUMP & RUN ############
dump.gsd("dihedral.gsd",
               period=100,
               group=group.all(),
               static=[],
               overwrite=True);
run(1e6);
