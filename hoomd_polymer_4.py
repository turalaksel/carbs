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

bead_positions = [(b.position[0][0],b.position[0][1],b.position[0][2]) for b in chains[0].nucleotides[0].beads]
bead_positions

com = np.average(np.asarray(bead_positions)[:,:3], axis=0)
bead_positions = np.asarray(bead_positions) - com

# Start HOOMD
context.initialize("");

snapshot = data.make_snapshot(N = n_nts,
                              box = hoomd.data.boxdim(Lx=120, Ly=120, Lz=120),
                              particle_types = types,
                              bond_types = ['polymer']);

snapshot.particles.position[:] = com_positions

snapshot.particles.typeid[:] = [0];

bonds = [[n, min(n+1, n_nts)] for n in range(0, n_nts - 1, 1)]

snapshot.bonds.resize(n_nts - 1);
snapshot.bonds.group[:] = bonds

system = init.read_snapshot(snapshot);

system.particles.types.add('A');

rigid = hoomd.md.constrain.rigid();


#can we get the quaternion of each base in the json?
rigid.set_param('C1',
                types=['A']*18,
                positions = bead_positions
                # positions=[(-0.4,0,0),(-0.3,0,0),(-0.2,0,0),(-0.1,0,0),
                #            (0.1,0,0),(0.2,0,0),(0.3,0,0),(0.4,0,0)]
                           );

rigid.create_bodies()

nl = hoomd.md.nlist.cell();

lj = hoomd.md.pair.lj(r_cut=2.0, nlist=nl)
lj.set_params(mode='shift')

lj.pair_coeff.set(['C1', 'A'], ['C1', 'A'], epsilon=1.0, sigma=1.0)


harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=10.0, r0=5);

hoomd.md.integrate.mode_standard(dt=0.005);

rigid = group.rigid_center();
hoomd.md.integrate.langevin(group=rigid, kT=0.5, seed=42);

hoomd.dump.gsd("trajectory_snap.gsd",
               period=1e3,
               group=hoomd.group.all(),
               overwrite=True);
hoomd.run(10e4);
