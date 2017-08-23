import numpy as np
import math
import quaternion as quat

###########################
# Useful Vector Functions #
###########################

def angle(p1,p2,p3,degrees=False):
  v1 = p2 - p1
  v2 = p3 - p2
  ang = np.arccos(np.dot(v2,v1)/(magvect(v1)*magvect(v2)))
  if degrees == True:
      ang = ang*180/np.pi
  return ang

def magvect(v):
  magnitude = np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return magnitude

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

def orthogonal(u):
    orth = np.random.randn(3)     # take a random vector
    orth -= np.dot(orth,u) * u    # make it orthogonal to k
    orth /= np.linalg.norm(orth)  # normalize it
    return(orth)

def calculateCoM(list_of_positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the center of mass (vector3) of the system of particles
    Assumes equal masses
    '''
    center = np.average(np.asarray(positions)[:,:3], axis=0)
    return(center)

def calculateMomentInertia(list_of_positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the moment of inertia (vector3) of the system of particles
    Assumes equal masses
    '''
    inertia = np.array([0., 0., 0.])
    center = calculateCoM(positions)
    for p, pos in enumerate(positions):
        shifted_pos = pos - center
        new_inertia = np.multiply(shifted_pos, shifted_pos)
        inertia = inertia + new_inertia
    #re-scale particle masses so that body is not hugely slow
    #this needs to be tested
    inertia /= p
    return(inertia)

###############################
# Useful Quaternion Functions #
###############################
def vectorQuaternion(u,v):
    cross = np.cross(u,v)
    if np.dot(u, v) < -0.999999:
        orthogonal = findOrthogonal(u)
        my_quaternion = quat.quaternion(0,*orthogonal)
    else:
        my_quaternion = quat.quaternion(0,*cross)
        my_quaternion.w = math.sqrt(np.linalg.norm(u)**2 * np.linalg.norm(v)**2) + np.dot(u,v)
    my_quaternion = quat.quaternion.normalized(my_quaternion)
    return(my_quaternion)

def testQuaternion(u,my_quat):
    rotated_vector = my_quat * quat.quaternion(*u) * np.conjugate(my_quat)
    return([rotated_vector.x,rotated_vector.y,rotated_vector.z])

#from https://stackoverflow.com/questions/16648452/calculating-quaternion-for-transformation-between-2-3d-cartesian-coordinate-syst
def systemQuaternion(lst1,lst2,matchlist=None):
    '''
    Given 2 lists of 3 orthogonal vectors,
    find the quaternion transformation bringing list2 into list1

    '''
    if not matchlist:
         matchlist=range(len(lst1))
    M=np.matrix([[0,0,0],[0,0,0],[0,0,0]])

    for i,coord1 in enumerate(lst1):
         x=np.matrix(np.outer(coord1,lst2[matchlist[i]]))
         M=M+x

    N11=float(M[0][:,0]+M[1][:,1]+M[2][:,2])
    N22=float(M[0][:,0]-M[1][:,1]-M[2][:,2])
    N33=float(-M[0][:,0]+M[1][:,1]-M[2][:,2])
    N44=float(-M[0][:,0]-M[1][:,1]+M[2][:,2])
    N12=float(M[1][:,2]-M[2][:,1])
    N13=float(M[2][:,0]-M[0][:,2])
    N14=float(M[0][:,1]-M[1][:,0])
    N21=float(N12)
    N23=float(M[0][:,1]+M[1][:,0])
    N24=float(M[2][:,0]+M[0][:,2])
    N31=float(N13)
    N32=float(N23)
    N34=float(M[1][:,2]+M[2][:,1])
    N41=float(N14)
    N42=float(N24)
    N43=float(N34)

    N=np.matrix([[N11,N12,N13,N14],\
                  [N21,N22,N23,N24],\
                  [N31,N32,N33,N34],\
                  [N41,N42,N43,N44]])

    values,vectors=np.linalg.eig(N)
    w=list(values)
    mw=max(w)
    my_quat= vectors[:,w.index(mw)]
    my_quat=np.array(my_quat).reshape(-1,).tolist()
    return quat.quaternion(*my_quat)
