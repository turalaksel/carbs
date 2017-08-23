# example dihedral test for hoomd v2.0.0
# Modified from Michael Howard's version
# see https://groups.google.com/d/msg/hoomd-users/CcpdtFg55M0/4jcGDPCyBAAJ

import numpy as np
from hoomd import *
from hoomd import md

context.initialize("");

# initialize 4 particles to the trans configuration with one diheral
s = data.make_snapshot(N=4, box=data.boxdim(L=50.0), improper_types=['I'])
s.particles.position[0] = [-1., 1., 0.]
s.particles.position[1] = [-1., 0., 0.]
s.particles.position[2] = [1., 0., 0.]
s.particles.position[3] = [1., -1., 0.]

s.impropers.resize(1)
s.impropers.group[0] = [0,1,2,3]

sys = init.read_snapshot(s)

# setup a standard harmonic diheral
improper_harm = md.improper.harmonic()
improper_harm.improper_coeff.set('I', k=50.0, chi=3.1415/2.)

# don't actually integrate anything, just testing the energies
md.integrate.mode_standard(dt=0.000)
all = group.all()
md.integrate.nve(group=all)

# sweep through dihedral angles and log the energy starting from the trans state
log = analyze.log(filename=None, quantities=['potential_energy'], period=1)
with open('improper_energy.log','w') as f:
    f.write('# phi U(phi)\n')
    for phi in np.linspace(-np.pi, np.pi, 360):
        p = sys.particles[3]
        p.position = (1., np.cos(3.1415-phi), np.sin(3.1415-phi))
        # p.position = (1. + np.cos(phi), np.sin(phi), 0)
        del p

        run(1)

        U = log.query('potential_energy')
        f.write('%f %f\n' % (phi, U))
