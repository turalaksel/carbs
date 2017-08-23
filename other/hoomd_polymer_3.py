import particlesFromPDB as fromPDB
from hoomd import *
import hoomd.md
import numpy as np

# Get chains of particles from PDB using particlesFromPDB.py
chains = fromPDB.getChainsFromPDB('doubleStrand.pdb')
n_nts = len(chains[0].nucleotides)
com_positions = [fromPDB.nucleotideCofM(n) for n in chains[0].nucleotides]
types = ['C1'.format(i) for i in range(0, n_nts)]

# shift positions to center of box
com = np.average(np.asarray(com_positions)[:,:3], axis=0)
com_positions = np.asarray(com_positions) - com

# Start HOOMD
context.initialize("");

uc = hoomd.lattice.unitcell(N = n_nts,
                            a1 = [120, 0,   0],
                            a2 = [0,    120, 0],
                            a3 = [0,    0,   120],
                            dimensions = 3,
                            position = com_positions,
                            type_name = types,
                            mass = [1.0]*n_nts,
                            moment_inertia = [[1,1,1]]*n_nts,
                            orientation = [[1, 0, 0, 0]]*n_nts
                            );

system = hoomd.init.create_lattice(unitcell=uc, n=[1,1,1]);

system.particles.types.add('A');

rigid = hoomd.md.constrain.rigid();
rigid.set_param('C1',
                types=['A']*8,
                positions=[(-0.4,0,0),(-0.3,0,0),(-0.2,0,0),(-0.1,0,0),
                           (0.1,0,0),(0.2,0,0),(0.3,0,0),(0.4,0,0)]);

rigid.create_bodies()

nl = hoomd.md.nlist.cell();

lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')

lj.pair_coeff.set(['C1', 'A'], ['C1', 'A'], epsilon=1.0, sigma=1.0)

hoomd.md.integrate.mode_standard(dt=0.005);

rigid = group.rigid_center();
hoomd.md.integrate.langevin(group=rigid, kT=1.0, seed=42);

hoomd.dump.gsd("trajectory_hard.gsd",
               period=1e3,
               group=hoomd.group.all(),
               overwrite=True);
hoomd.run(1e4);
