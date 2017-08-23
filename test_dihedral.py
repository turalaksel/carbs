# example dihedral test for hoomd v2.0.0
# Modified from Michael Howard's version
# see https://groups.google.com/d/msg/hoomd-users/CcpdtFg55M0/4jcGDPCyBAAJ

import numpy as np
from hoomd import *
from hoomd import md

context.initialize("");

# initialize 4 particles to the trans configuration with one diheral
s = data.make_snapshot(N=4, box=data.boxdim(L=50.0), dihedral_types=['D'], angle_types=['A'])
s.particles.position[0] = [-1., 1., 0.]
s.particles.position[1] = [-1., 0., 0.]
s.particles.position[2] = [1., 0., 0.]
s.particles.position[3] = [1., -1., 0.]

s.dihedrals.resize(1)
s.dihedrals.group[0] = [0,1,2,3]

s.angles.resize(1)
s.angles.group[0] = [1,2,3]

sys = init.read_snapshot(s)

# setup a standard harmonic diheral
def harmonicAngle(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)

dtable = md.dihedral.table(width=1000)
dtable.dihedral_coeff.set('D', func=harmonicAngle, coeff=dict(kappa=50., theta0=-3.1415/2.))

dharm = md.dihedral.harmonic()
dharm.dihedral_coeff.set('D', k=50.0, d=-1, n=1)

# aharm = md.angle.harmonic()
# aharm.angle_coeff.set('A', k=50.0, t0=0)

# don't actually integrate anything, just testing the energies
md.integrate.mode_standard(dt=0.000)
all = group.all()
md.integrate.nve(group=all)

# sweep through dihedral angles and log the energy starting from the trans state
log = analyze.log(filename=None, quantities=['potential_energy'], period=1)
with open('dihedral_energy.log','w') as f:
    f.write('# phi U(phi)\n')
    for phi in np.linspace(-np.pi, np.pi, 360):
        p = sys.particles[3]
        # p.position = (1., np.cos(phi), np.sin(phi))
        p.position = (1., np.cos(phi), np.sin(phi))
        del p

        run(1)

        U = log.query('potential_energy')
        f.write('%f %f\n' % (phi, U))
