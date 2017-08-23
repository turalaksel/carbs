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
for line in open('doubleStrand.pdb'):
    line = line.split()
    particle_type = line[0]
    chain = 0
    if particle_type == 'ATOM':
        particle_number = int(line[1])
        types.append(int(line[5]) - 1)
        positions.append(line[6:9])
    elif particle_type == 'TEM':
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
