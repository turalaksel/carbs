import numpy as np
import math

###########################
# Useful Vector Functions #
###########################

def calculate_comass(list_of_positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the center of mass (vector3) of the system of particles
    Assumes equal masses
    '''
    center = np.average(np.asarray(list_of_positions)[:,:3], axis=0)
    return(center)

def calculate_moment_inertia(list_of_positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the moment of inertia (vector3) of the system of particles
    Assumes equal masses
    '''
    inertia = np.array([0., 0., 0.])
    center = calculate_comass(list_of_positions)
    for p, pos in enumerate(list_of_positions):
        shifted_pos = pos - center
        new_inertia = np.multiply(shifted_pos, shifted_pos)
        inertia = inertia + new_inertia
    
    return(inertia)

def eular_distance(particle_coor_1,particle_coor_2):
    '''
    Measure the eular distance between two coordinates
    '''
    diff_vec = np.array(particle_coor_2) - np.array(particle_coor_1)
    
    return np.sqrt(np.sum(diff_vec**2))
