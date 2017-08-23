import numpy as np
import particlesFromPDB as fromPDB

chains = fromPDB.getChainsFromPDB("/Users/damasceno/Desktop/2017-08-05/oxDNA.pdb")

# angle stuff
def magvect(v):
  magnitude = np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return magnitude

# angle between vectors
def angle(p1,p2,p3,degrees=False):
  v1 = p2 - p1
  v2 = p3 - p2
  ang = np.arccos(np.dot(v2,v1)/(magvect(v1)*magvect(v2)))
  if degrees == True:
      ang = ang*180/np.pi
  return ang

def projection(v1,v2):
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    v1_hat = v1/v1_norm
    v2_hat = v2/v2_norm
    cosTheta = np.dot(v1,v2)/(v1_norm * v2_norm)
    projection = v1_norm * cosTheta * v2_hat
    return(projection)

def rejection(v1,v2):
    rejection = v1 - projection(v1,v2)
    return(rejection)


c = []
b = []
a = []
for i in range(0, len(chains[0].nucleotides)):
    nucleotide = chains[0].nucleotides[i]
    c.append(nucleotide.beads[0].position[0])
    b.append(nucleotide.beads[1].position[0])

c = np.asarray(c)
b = np.asarray(b)

for i in range(len(c)):
    if i < len(c) - 1:
        b1_ref = b[i] - c[i]
        c2_ref = c[i+1] - c[i]
        new_c_ref = rejection(c2_ref, b1_ref)
        aa = np.cross(new_c_ref,b1_ref)
    else:
        b1_ref = b[i] - c[i]
        c2_ref = c[i-1] - c[i]
        new_c_ref = rejection(c2_ref, b1_ref)
        aa = np.cross(new_c_ref,b1_ref)
    a.append(list(aa))

a = np.asarray(a)
a = a + c

for i in range(len(c)-2):
    angle1 = angle(c[i], c[i+1], a[i+1],False)
    angle2 = angle(c[i], c[i+1], b[i+1],False)
    print(angle1,angle2)

def dihedral(p1,p2,p3,p4,degrees=False):
    v1 = p1 - p2
    v2 = p3 - p2
    v3 = p4 - p3
    v2m = -v2
    a = np.cross(v1, v2m)
    b = np.cross(v3, v2m)

    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)
    v2m_norm = np.linalg.norm(v2m)

    cosPhi = np.dot(a, b) / (a_norm * b_norm)
    sinPhi = v2m_norm * np.dot(a, v3) / (a_norm * b_norm)

    dihedral_angle = np.arctan(sinPhi/ cosPhi)
    return(dihedral_angle)

# for i in range(len(c)-1):
#     dihedral1 = dihedral(b[i],c[i],c[i+1],b[i+1],False)
#     dihedral2 = dihedral(b[i],c[i],a[i],c[i+1],False)
#     dihedral3 = dihedral(a[i],c[i],b[i],c[i+1],False)
#     print(dihedral1, dihedral2, dihedral3)

for i in range(len(c)-2):
    angle1 = dihedral(c[i], c[i+1], b[i+1], c[i+2],False)
    angle2 = dihedral(c[i], c[i+1], a[i+1], c[i+2],False)
    print(angle1, angle2)

len(c)
