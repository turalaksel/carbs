from hoomd import *
from hoomd import md
import numpy as np
import random
import particlesFromPDB as fromPDB

# Start HOOMD
context.initialize("");

num_part = 20 #even !!

snapshot = data.make_snapshot(N = num_part,
                              box = data.boxdim(Lx=30, Ly=30, Lz=num_part + 10),
                              particle_types=['C','a','b'], #center, sidea, sideb
                              bond_types = ['polymer'],
                              angle_types = [angle1, angle2],
                              dihedral_types = ['dihedral','dihedral2','dihedral3'],
                              improper_types = []);

#positions
def ran():
    return(round(random.uniform(-0.1, 0.1),3))
part_pos = [[0.,0.,i] for i in range(-int(num_part/2),int(num_part/2))]
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
rigid.set_param('C', types=['A','B'], positions = [[1.0,0.,0.],[0,1.0,0.]]);
rigid.create_bodies()

for i in range(num_part):
    system.particles[0].diameter = 1.0
    system.particles[num_part + 2*i].diameter = 0.5
    system.particles[num_part + 2*i + 1].diameter = 0.1

#fene / harmonic
harmonic1 = md.bond.harmonic()
harmonic1.bond_coeff.set('polymer', k=50.0, r0=1.2)

#angle for backbone
for i in range(num_part-1):
    system.angles.add('angle1', i, i+1, i+2)

harmonic2 = md.angle.harmonic()
harmonic2.angle_coeff.set('angle', k=300.0, t0=2.6)


# #dihedral
# def harmonicAngle(theta, kappa, theta0):
#    V = 0.5 * kappa * (theta-theta0)**2;
#    F = -kappa*(theta-theta0);
#    return (V, F)
# for i in range(num_part-1):
#     system.dihedrals.add('dihedral', i, num_part + 2*i, i+1, num_part + 2*(i+1))
# for i in range(num_part-1): #
#     system.dihedrals.add('dihedral2', num_part + 2*i, i+1, i, num_part + 2*i+1)
# for i in range(0,num_part-1):
#     system.dihedrals.add('dihedral3', num_part + 2*i + 1, i+1, i, num_part + 2*i)
#
# dtable = md.dihedral.table(width=1000)
# dtable.dihedral_coeff.set('dihedral', func=harmonicAngle, coeff=dict(kappa=10, theta0=0.))
# dtable.dihedral_coeff.set('dihedral2', func=harmonicAngle, coeff=dict(kappa=50, theta0=-1.57))
# dtable.dihedral_coeff.set('dihedral3', func=harmonicAngle, coeff=dict(kappa=0, theta0=1.57))

# LJ interactions
wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
wca.set_params(mode='shift')
wca.pair_coeff.set('C', 'C', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))
wca.pair_coeff.set('C', 'A', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))
wca.pair_coeff.set('C', 'B', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))
wca.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))
wca.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))
wca.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))

md.integrate.mode_standard(dt=0.005, aniso=True);

rigid = group.rigid_center();
md.integrate.langevin(group=rigid, kT=0.1, seed=42);

dump.gsd("dihedral.gsd",
               period=10,
               group=group.all(),
               overwrite=True);
run(100000);
