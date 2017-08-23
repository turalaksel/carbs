import particlesFromPDB as fromPDB
from hoomd import *
import hoomd.md

chains = fromPDB.getChainsFromPDB('doubleStrand.pdb')

chains[1].nucleotides[0].beads[0].position

context.initialize("");

uc = lattice.unitcell(N = len(pos_array),
                            a1 = [200., 0,   0],
                            a2 = [0,    200., 0],
                            a3 = [0,    0,   200.],
                            dimensions = 3,
                            position = pos_array,
                            type_name = ['A']*len(pos_array),
                            mass = [1.]*len(pos_array),
                            moment_inertia = [[1.,1.,1.]]*len(pos_array),
                            orientation = [[1, 0, 0, 0]]*len(pos_array) );

system = init.create_lattice(unitcell=uc, n=[1,1,1]);

nl = hoomd.md.nlist.cell()

lj = hoomd.md.pair.lj(r_cut=2**(1/6), nlist=nl)
lj.set_params(mode='shift')

lj.pair_coeff.set(['A', 'A'], ['A', 'A'], epsilon=0.0, sigma=1.0)

all = group.all()
hoomd.md.integrate.mode_standard(dt=0.002) #################### THEY USE 0.002
hoomd.md.integrate.nvt(group=all, tau=1.0, kT=0.05)

hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy',
                              'translational_kinetic_energy',
                              'rotational_kinetic_energy'],
                  period=100,
                  overwrite=True);

hoomd.dump.gsd("trajectory.gsd",
               period=100,
               group=hoomd.group.all(),
               overwrite=True);

hoomd.run(401)
