# example dihedral test for hoomd v2.0.0
# Modified from Michael Howard's version
# see https://groups.google.com/d/msg/hoomd-users/CcpdtFg55M0/4jcGDPCyBAAJ

import numpy as np
from hoomd import *
from hoomd import md

context.initialize("");

# initialize 4 particles to the trans configuration with one diheral
s = data.make_snapshot(N=3, box=data.boxdim(L=50.0), dihedral_types=['D'], angle_types=['A'])
s.particles.position[0] = [0., 0., 0.]
s.particles.position[1] = [1., 0., 0.]
s.particles.position[2] = [2., 0., 0.]

s.angles.resize(1)
s.angles.group[0] = [0,1,2]

sys = init.read_snapshot(s)

# setup a standard harmonic angle

aharm = md.angle.harmonic()
aharm.angle_coeff.set('A', k=50.0, t0=3*3.1415/4)

# don't actually integrate anything, just testing the energies
md.integrate.mode_standard(dt=0.000)
all = group.all()
md.integrate.nve(group=all)

# sweep through dihedral angles and log the energy starting from the trans state
log = analyze.log(filename=None, quantities=['potential_energy'], period=1)
energies=[]
positions=[]
with open('dihedral_energy.log','w') as f:
    f.write('# phi U(phi)\n')
    for phi in np.linspace(-2*np.pi, 2*np.pi, 360):
        p = sys.particles[2]
        # p.position = (1., np.cos(phi), np.sin(phi))
        p.position = (1. + np.cos(phi), np.sin(phi), 0)
        positions.append(p.position)
        del p

        run(1)

        U = log.query('potential_energy')
        energies.append(U)
        f.write('%f %f\n' % (phi, U))

import matplotlib.pyplot as plt

x = np.linspace(-np.pi, np.pi, 360);
y = np.asarray(energies)
plt.plot(x, y)
plt.show()
