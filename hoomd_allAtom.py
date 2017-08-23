'''
This program generates a HOOMD all-atom simulation
from a PDB file from COSM
'''

import hoomd
hoomd.context.initialize("");

uc = hoomd.lattice.unitcell(N = 1,
                            a1 = [300., 0,   0],
                            a2 = [0,    300., 0],
                            a3 = [0,    0,   300.],
                            dimensions = 3,
                            position = positions_array,
                            type_name = types_array,
                            mass = masses_array,
                            # body = bodies_array,
                            moment_inertia = [[0,
                                               1/12*1.0*8**2,
                                               1/12*1.0*8**2]],
                            orientation = [[1, 0, 0, 0]]);

system = hoomd.init.create_lattice(unitcell=uc, n=[2,18,18]);

# first create the center of mass particles

# then attach the nucleotides to each CoM
