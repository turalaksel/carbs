import numpy as np

import argparse
from os import path
import sys
import cadnano
from cadnano.document import Document
import vector_tools as vTools

# Read nucleotides data from cadnano

app = cadnano.app()
doc = app.document = Document()
doc.readFile('input/6hb.json');

part = doc.activePart()
oligos = part.oligos()

#get all oligos (including scaffold)
oligos_sorted_by_length = sorted(oligos, key=lambda x: x.length(), reverse=True)
longest_oligo = oligos_sorted_by_length[0]
staple_oligos = oligos_sorted_by_length[1:]
oligos_array = [longest_oligo] + [staple for staple in staple_oligos]

##################################
#some helper functions for cadnano
##################################
def oligoHelperList(oligo):
    '''
    Given an oligo, returns the following list of useful properties
    (strand, strand direction, vh id, index)
    for each nucleotide. List is oriented 5' to 3' and ends back
    at 1st particle if oligo is circular.

    Args:
        oligo

    Returns:
        list containing strand, strand direction, vh id, and index

    '''

    generator = oligo.strand5p().generator3pStrand()
    oligo_helper_list = []
    for strand in generator:
        vh = strand.idNum()
        index_5p = strand.idx5Prime()
        index_3p = strand.idx3Prime()

        direction = (-1 + 2*strand.isForward()) #-1 if backwards, 1 if fwd
        for i in range(index_5p, index_3p + direction, direction):
            oligo_helper_list.append([strand, direction, vh, i])
    if oligo.isCircular() == True: #connect 1st and last atoms
        oligo_helper_list.append(oligo_helper_list[0])
    return(oligo_helper_list)

def findCoordinates(vh, index):
    '''
    Given a vh and a index, returns (x,y,z)
    for the backbone pts and sidechains fwd and rev
    '''
    axis_pts = part.getCoordinates(vh)[0][index]
    fwd_pts = part.getCoordinates(vh)[1][index]
    rev_pts = part.getCoordinates(vh)[2][index]
    return [fwd_pts, axis_pts, rev_pts]

class nucleotide:
    '''
    Fixed attributes of a nucleotide
    Attributes: index, position, strand, vh
    '''
    def __init__(self):
        self.direction  = 1 #1 is fwd, -1 is reverse
        self.index      = []
        self.position   = []
        self.quaternion = []
        self.strand	    = []
        self.vectors    = []
        self.vh  	    = []
    def add_index(self, index):
        self.index.append(index)
    def add_position(self, position):
        self.position.append(position)
    def add_quaternion(self, quat):
        self.quaternion.append(quat)
    def add_strand(self, strand):
        self.strand.append(strand)
    def add_vectors(self, vectors):
        self.vectors.append(vectors)
    def add_vh(self, vh):
        self.vh.append(vh)

def populateNucleotides(oligo):
    nucleotides_list = []
    for [strand, direction, vh, index] in oligoHelperList(oligo):
        coordinates = findCoordinates(vh, index)

        nucl = nucleotide()
        nucl.direction = direction
        nucl.add_index(index)
        nucl.add_position(coordinates[1]) #axis / side-chain position
        nucl.add_position(coordinates[1 + direction]) #fwd or rev backbone vector
        nucl.add_strand(strand)
        nucl.add_vh(vh)
        nucleotides_list.append(nucl)
    return(nucleotides_list)

def listOflistsOfNucleotides(oligos_array):
    nucl_list_list = []
    for oligo in oligos_array:
        nucleotides_list = populateNucleotides(oligo)
        nucl_list_list.append(nucleotides_list)
    return(nucl_list_list)

def generateVectorsandQuaternions(oligos_array):
    '''
    Given an array of oligos, fills in nucleotide attributes
    '''
    nucl_list_list = listOflistsOfNucleotides(oligos_array)

    for o, oligo in enumerate(nucl_list_list):
        [axis_0, backbone_0] = nucl_list_list[0][0].position
        [axis_1, backbone_1] = nucl_list_list[0][1].position

        base_vector_0 = axis_0 - backbone_0
        backbone_vector_0 = backbone_1 - backbone_0
        aux_vector_a_0 = np.cross(base_vector_0, backbone_vector_0)
        aux_vector_b_0 = np.cross(aux_vector_a_0, base_vector_0)
        # return 3 orthogonal vectors in nucleotide, for quaternion
        vect_list_0 = (base_vector_0/np.linalg.norm(base_vector_0), \
                     aux_vector_a_0/np.linalg.norm(aux_vector_a_0), \
                     aux_vector_b_0/np.linalg.norm(aux_vector_b_0))

        for n, nucl in enumerate(nucl_list_list[o]):
            [axis_1, backbone_1] = nucl_list_list[o][n].position
            if n < len(nucl_list_list[o]) - 1:
                [axis_2, backbone_2] = nucl_list_list[o][n + 1].position
                base_vector_1 = axis_1 - backbone_1
                backbone_vector_1 = backbone_2 - backbone_1
            elif n == len(nucl_list_list[o]) - 1:
                [axis_2, backbone_2] = nucl_list_list[o][n - 1].position
                base_vector_1 = axis_1 - backbone_1
                backbone_vector_1 = - (backbone_2 - backbone_1)

            aux_vector_a_1 = np.cross(base_vector_1, backbone_vector_1)
            aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
            vect_list_1 = (base_vector_1+np.array([0.00001,0,0])/np.linalg.norm(base_vector_1+np.array([0.00001,0,0])), \
                         aux_vector_a_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_a_1+np.array([0.00001,0,0])), \
                         aux_vector_b_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_b_1+np.array([0.00001,0,0])))

            nucl_quaternion = vTools.systemQuaternion(vect_list_0, vect_list_1)
            nucl.add_vectors(vect_list_1)
            nucl.add_quaternion(nucl_quaternion.w)
            nucl.add_quaternion(nucl_quaternion.x)
            nucl.add_quaternion(nucl_quaternion.y)
            nucl.add_quaternion(nucl_quaternion.z)
    return(nucl_list_list)

##################################
#start HOOMD code
##################################
from hoomd import *
from hoomd import md

# Start HOOMD
context.initialize("");

# calculate list of oligos, each as a list of nucleotides in 5'-3' direction
list_of_list_of_nucleotides = generateVectorsandQuaternions(oligos_array)

# calculate positions and save them into a flattened Nx3 array
nucl_positions = [[[nucl.position[1]] for nucl in nucl_chains] for nucl_chains in list_of_list_of_nucleotides]
nucl_positions = np.concatenate(np.concatenate(nucl_positions))

nucl_quaternions = [[nucl.quaternion for nucl in nucl_chains] for nucl_chains in list_of_list_of_nucleotides]
nucl_quaternions = np.concatenate(nucl_quaternions)

num_oligos = len(list_of_list_of_nucleotides)
total_num_nucl = len(nucl_positions)

snapshot = data.make_snapshot(N = total_num_nucl,
                              box = data.boxdim(Lx=30, Ly=30, Lz=30),
                              particle_types=['backbone','sidechain','aux'],
                              bond_types = ['backbone'],
                              dihedral_types = ['dihedral1', \
                                                'dihedral21',\
                                                'dihedral22',\
                                                'dihedral31',\
                                                'dihedral32']);

# particle positions, types and moments of inertia
snapshot.particles.position[:] = nucl_positions
snapshot.particles.moment_inertia[:] = [[1.,1.,1.]] #not correct. fix it
snapshot.particles.typeid[:] = [0];

# Backbone bonds
bonds = []
i = 0
for chain in list_of_list_of_nucleotides:
    for n in range(i, i + len(chain) - 1):
        bonds.append([n, n+1])
    i += len(chain)

snapshot.bonds.resize(total_num_nucl - num_oligos)
snapshot.bonds.group[:] = bonds

# Read the snapshot and create neighbor list
system = init.read_snapshot(snapshot);
nl = md.nlist.cell();

############ BONDS ############
#rigid
nucl0 = list_of_list_of_nucleotides[0][0]
rigid = md.constrain.rigid();
rigid.set_param('backbone', types=['sidechain','aux'], positions = [0.9*nucl0.vectors[0][0], 0.4*nucl0.vectors[0][1]]);
rigid.create_bodies()

#fene / harmonic
harmonic1 = md.bond.harmonic()
harmonic1.bond_coeff.set('backbone', k=10.0, r0=0.75)

# dihedrals
def harmonicAngle(theta, kappa, theta0):
   V = 0.5 * kappa * (theta-theta0)**2;
   F = -kappa*(theta-theta0);
   return (V, F)

# for chain in range(len(list_of_list_of_nucleotides)):
#     for nucl in range(len(list_of_list_of_nucleotides[chain]) - 2):
#         chain_length = len(list_of_list_of_nucleotides[chain])
#         sidechain_1 = total_num_nucl + 2*chain*chain_length + 2*nucl
#         aux_1 = total_num_nucl + 2*chain*chain_length + 2*nucl + 1
#         bckbone_1 = chain*chain_length + nucl
#         sidechain_2 = total_num_nucl + 2*chain*chain_length + 2*nucl + 2
#         aux_2 = total_num_nucl + 2*chain*chain_length + 2*nucl + 3
#         bckbone_2 = chain*chain_length + nucl + 1
#
#         system.dihedrals.add('dihedral1',  aux_1, bckbone_1, bckbone_2, aux_2)
#         system.dihedrals.add('dihedral21', aux_1, bckbone_1, sidechain_1, bckbone_2)
#         system.dihedrals.add('dihedral22', aux_2, bckbone_2, sidechain_2, bckbone_1)
#         system.dihedrals.add('dihedral31', sidechain_1, bckbone_1, aux_1, bckbone_2)
#         system.dihedrals.add('dihedral32', sidechain_2, bckbone_2, aux_2, bckbone_1)
#
# dtable = md.dihedral.table(width=1000)
# dtable.dihedral_coeff.set('dihedral1',  func=harmonicAngle, coeff=dict(kappa=20, theta0=-0.33))
# dtable.dihedral_coeff.set('dihedral21', func=harmonicAngle, coeff=dict(kappa=20, theta0=+1.23)) #1.33
# dtable.dihedral_coeff.set('dihedral22', func=harmonicAngle, coeff=dict(kappa=20, theta0=-1.23))
# dtable.dihedral_coeff.set('dihedral31', func=harmonicAngle, coeff=dict(kappa=20, theta0=-1.67)) #1.57
# dtable.dihedral_coeff.set('dihedral32', func=harmonicAngle, coeff=dict(kappa=20, theta0=+1.67))


# fix particle quaternions
p=0
for chain in list_of_list_of_nucleotides:
    for nucl in chain:
        nucl_quaternion = nucl.quaternion
        system.particles[p].orientation = nucl_quaternion
        p += 1

# fix diameters for vizualization
for i in range(0, total_num_nucl):
    system.particles[i].diameter = 0.8
for i in range(total_num_nucl, len(system.particles), 2):
    system.particles[i].diameter = 0.3
    system.particles[i + 1].diameter = 0.1


########## INTERACTIONS ############
# LJ interactions
wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
wca.set_params(mode='shift')
wca.pair_coeff.set('backbone', 'backbone',   epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
wca.pair_coeff.set('backbone', 'sidechain',  epsilon=0.0, sigma=0.665, r_cut=0.665*2**(1/6))
wca.pair_coeff.set('sidechain', 'sidechain', epsilon=0.0, sigma=0.330, r_cut=0.330*2**(1/6))
wca.pair_coeff.set('backbone', 'aux',        epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
wca.pair_coeff.set('aux', 'sidechain',       epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
wca.pair_coeff.set('aux', 'aux',             epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))

########## INTEGRATION ############
md.integrate.mode_standard(dt=0.003, aniso=True);
rigid = group.rigid_center();
md.integrate.langevin(group=rigid, kT=0.2, seed=42);
########## DUMP & RUN ############
dump.gsd("output/6hb.gsd",
               period=1,
               group=group.all(),
               static=[],
               overwrite=True);
run(1000);









#
#
