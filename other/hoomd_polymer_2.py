import particlesFromPDB as fromPDB
from hoomd import *
import hoomd.md
import numpy as np

# Get chains of particles from PDB using particlesFromPDB.py
chains = fromPDB.getChainsFromPDB('doubleStrand.pdb')
n_nts = len(chains[0].nucleotides)
com_positions = [fromPDB.nucleotideCofM(n) for n in chains[0].nucleotides]
types = ['C{}'.format(i) for i in range(0, n_nts)]

# shift positions to center of box
com = np.average(np.asarray(com_positions)[:,:3], axis=0)
com_positions = np.asarray(com_positions) - com

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

snapshot.particles.velocity[:] = np.random.normal(0.0,
    np.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

system = init.read_snapshot(snapshot);

system.particles.types.add('A');

rigid = hoomd.md.constrain.rigid();
rigid.set_param('C1',
                types=['A']*8,
                positions=[(-0.4,0,0),(-0.3,0,0),(-0.2,0,0),(-0.1,0,0),
                           (0.1,0,0),(0.2,0,0),(0.3,0,0),(0.4,0,0)]);

rigid.create_bodies()


nl = hoomd.md.nlist.cell();

dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.8, seed=1);

dpd.pair_coeff.set('A', 'A', A=25.0, gamma = 1.0);

nl.reset_exclusions(exclusions = []);

harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=100.0, r0=0);

hoomd.md.integrate.mode_standard(dt=0.01);

all = hoomd.group.all();
hoomd.md.integrate.nve(group=all);

hoomd.dump.gsd("trajectory.gsd", period=10, group=all, overwrite=True);

hoomd.run(5e4);
