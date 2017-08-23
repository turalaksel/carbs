'''
This program transforms a PDB file generated from COSM
into a GSD file that can be read by HOOMD
'''

import gsd
import gsd.hoomd
import gsd.pygsd
import numpy as np

# Start the GSD file
t = gsd.hoomd.open(name='test.gsd', mode='wb')

# Define the variables to append
positions = []
types = []

#read the pdb
previous_body = 1
group_particle_number = []
group_positions = []
group_types = []
cm = []
for line in open('doubleStrand.pdb'):
    line = line.split()
    is_atom = line[0]
    chain = 0

    if is_atom == 'ATOM':
        this_body = int(line[5]) - 1
        if this_body == previous_body:
            # group_particle_number.append(int(line[1]))
            group_types.append(int(line[5]) - 1)
            group_positions.append(line[6:9])
        else: #we found another body. Append all group particles
            # particle_number.append(group_particle_number)
            types.append(group_types)
            positions_array = np.asarray(group_positions, dtype=np.float32)
            cm.append(np.average(positions_array[:,:3], axis=0))

            group_particle_number = []
            group_positions = []
            group_types = []

    elif is_atom == 'TEM':
        chain += 1

# correct for center of mass
pos_array = np.asarray(positions, dtype=np.float32)
CM = np.average(pos_array[:,:3], axis=0)
pos_array = pos_array - CM

s = gsd.hoomd.Snapshot()
s.configuration.step = 0
s.configuration.box = [100, 100, 100, 0, 0, 0]

s.particles.N = len(pos_array)
s.particles.position = pos_array

t.append( s )
