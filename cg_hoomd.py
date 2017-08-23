from hoomd import *
from hoomd import md
import numpy as np

# Start HOOMD
context.initialize("");

bps_per_polymer = 60 #even !!
num_poly = 2
num_part = bps_per_polymer * num_poly

snapshot = data.make_snapshot(N = bps_per_polymer,
                              box = data.boxdim(Lx=bps_per_polymer/2, Ly=bps_per_polymer, Lz=bps_per_polymer + 30),
                              particle_types=['C','a','b'], #center, sidea, sideb
                              bond_types = ['polymer'],
                              dihedral_types = ['dihedral1','dihedral21','dihedral22','dihedral31','dihedral32'],
                              improper_types = []);

############ INIT ############
#positions
part_pos = [[0.,0.,i] for i in range(-int(bps_per_polymer/2),int(bps_per_polymer/2))]
snapshot.particles.position[:] = part_pos
#inertia
snapshot.particles.moment_inertia[:] = [[1./3.,1./3.,1./3.]]*bps_per_polymer
#type
snapshot.particles.typeid[:] = [0];
#bonds
bonds = [[n, min(n+1, bps_per_polymer)] for n in range(0, bps_per_polymer - 1, 1)]
snapshot.bonds.resize(bps_per_polymer - 1);
snapshot.bonds.group[:] = bonds
snapshot.replicate(2,1,1);
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
harmonic1.bond_coeff.set('polymer', k=50.0, r0=0.7)

#dihedrals
def harmonicAngle(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)

for n in range(num_poly):
    for i in range(bps_per_polymer-1):
        a1 = num_part + 2*n*bps_per_polymer + 2*i
        b1 = num_part + 2*n*bps_per_polymer + 2*i + 1
        c1 = n*bps_per_polymer + i
        a2 = num_part + 2*n*bps_per_polymer + 2*i + 2
        b2 = num_part + 2*n*bps_per_polymer + 2*i + 3
        c2 = n*bps_per_polymer + i + 1
        system.dihedrals.add('dihedral1', b1, c1, c2, b2)
        system.dihedrals.add('dihedral21', b1, c1, a1, c2)
        system.dihedrals.add('dihedral22', b2, c2, a2, c1)
        system.dihedrals.add('dihedral31', a1, c1, b1, c2)
        system.dihedrals.add('dihedral32', a2, c2, b2, c1)

dtable = md.dihedral.table(width=1000)
dtable.dihedral_coeff.set('dihedral1',  func=harmonicAngle, coeff=dict(kappa=20, theta0=-0.33))
dtable.dihedral_coeff.set('dihedral21', func=harmonicAngle, coeff=dict(kappa=20, theta0=+1.23)) #1.33
dtable.dihedral_coeff.set('dihedral22', func=harmonicAngle, coeff=dict(kappa=20, theta0=-1.23))
dtable.dihedral_coeff.set('dihedral31', func=harmonicAngle, coeff=dict(kappa=20, theta0=-1.67)) #1.57
dtable.dihedral_coeff.set('dihedral32', func=harmonicAngle, coeff=dict(kappa=20, theta0=+1.67))

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
