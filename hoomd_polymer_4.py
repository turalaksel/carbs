import particlesFromPDB as fromPDB
from hoomd import *
import hoomd.md
import numpy as np

# Get chains of particles from PDB using particlesFromPDB.py
chains = fromPDB.getChainsFromPDB('double_strand.pdb')
n_nts = len(chains[0].nucleotides)
types = ['C1'.format(i) for i in range(0, n_nts)]

# shift positions to center of box
com_positions = [fromPDB.nucleotideCofM(n) for n in chains[0].nucleotides]
com = np.average(np.asarray(com_positions)[:,:3], axis=0)
com_positions = np.asarray(com_positions) - com

bead_positions = [(b.position[0][0],b.position[0][1],b.position[0][2]) \
    for b in chains[0].nucleotides[1].beads]

com = np.average(np.asarray(bead_positions)[:,:3], axis=0)
bead_positions = np.asarray(bead_positions) - com

# Start HOOMD
context.initialize("");

snapshot = data.make_snapshot(N = n_nts,
                              box = hoomd.data.boxdim(Lx=150, Ly=150, Lz=150),
                              particle_types = types,
                              bond_types = ['polymer']);

snapshot.particles.position[:] = com_positions

snapshot.particles.moment_inertia[:] = [[100,100,100]]*n_nts

snapshot.particles.typeid[:] = [0];

bonds = [[n, min(n+1, n_nts)] for n in range(0, n_nts - 1, 1)]

snapshot.bonds.resize(n_nts - 1);
snapshot.bonds.group[:] = bonds

system = init.read_snapshot(snapshot);

system.particles.types.add('A');

rigid = hoomd.md.constrain.rigid();


#can we get the quaternion of each base in the json?
rigid.set_param('C1',
                types=['A']*22,
                positions = bead_positions
                           );

rigid.create_bodies()

nl = hoomd.md.nlist.cell();

lj = hoomd.md.pair.lj(r_cut=2.0, nlist=nl)
lj.set_params(mode='shift')

lj.pair_coeff.set(['C1', 'A'], ['C1', 'A'], epsilon=1.0, sigma=1.0)

nl.reset_exclusions(exclusions = []);

harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=20.0, r0=7);

def harmonic(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)

dtable = hoomd.md.dihedral.table(width=1000)
dtable.dihedral_coeff.set('polymer', func=harmonic, coeff=dict(kappa=3000, theta_0=0.))


hoomd.md.integrate.mode_standard(dt=0.005, aniso=True);

rigid = group.rigid_center();
hoomd.md.integrate.langevin(group=rigid, kT=0.5, seed=42);

hoomd.dump.gsd("trajectory_snap.gsd",
               period=1e3,
               group=hoomd.group.all(),
               overwrite=True);
hoomd.run(10e4);
