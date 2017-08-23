import numpy as np
import argparse
import sys
import cadnano
from os import path
from cadnano.document import Document
import vector_tools as vTools #needed for quaternion math and other vector calculations

####### USER DEFINED VARIABLES########
INPUT_FILENAME = 'input/tripod.json'
RELAX = True
######################################

##################################
#  Helper functions for cadnano  #
##################################
class Nucleotide:
    '''
    Fixed attributes of a nucleotide
    Attributes: index, position, strand, vh
    '''
    def __init__(self):
        self.direction  = 1  #1 is fwd, -1 is reverse
        self.global_pts = None # backbone, sidechain and aux points in global frame
        self.index      = None # z position in cadnano's unit
        self.body       = None # body this nucl belongs to (for relaxation)
        self.position   = None # positions of axis (sidechain) and backbone ptcls
        self.quaternion = None
        self.strand	    = None # strand #
        self.vectors    = None # orth vectors in the body ref. frame for quat calc
        self.vh  	    = None

class Body:
    '''
    Fixed attributes of a body.
    A body is a combination of neighboring vhs, making up a
    collection of nucleotides that move together
    Attributes: vhs, com_position, com_quaternion, nucleotides
    '''
    def __init__(self):
        self.com_position   = None
        self.com_quaternion = None
        self.moment_inertia = None
        self.nucleotides    = []   # nucleotides composing the body
        self.vhs            = set() # vhs composing the body

    def add_nucleotide(self, nucl):
        self.nucleotides.append(nucl)
    def add_vh(self, vh):
        self.vhs.add(vh)

def oligoHelperList(oligo):
    '''
    Given an oligo, returns the following list of useful properties
    (strand, strand direction, vh id, index)
    for each nucleotide. List is oriented 5' to 3' and ends back
    at 1st particle if oligo is circular.
    '''

    generator = oligo.strand5p().generator3pStrand()
    oligo_helper_list = []
    for strand in generator:
        strand_helper_list = []
        vh = strand.idNum()
        index_5p = strand.idx5Prime()
        index_3p = strand.idx3Prime()

        direction = (-1 + 2*strand.isForward()) #-1 if backwards, 1 if fwd
        for i in range(index_5p, index_3p + direction, direction):
            strand_helper_list.append([strand, direction, vh, i])
        oligo_helper_list.append(strand_helper_list)
    return(oligo_helper_list)

def findCoordinates(vh, index):
    '''
    Given a vh and a index, returns (x,y,z)
    for the sidechain pts and backbones fwd and rev
    '''
    axis_pts = part.getCoordinates(vh)[0][index]
    fwd_pts = part.getCoordinates(vh)[1][index]
    rev_pts = part.getCoordinates(vh)[2][index]
    return [fwd_pts, axis_pts, rev_pts]

def populateNucleotides(oligo):
    strand_list = []
    for strands in oligoHelperList(oligo):
        nucleotides_list = []
        for [strand, direction, vh, index] in strands:
            coordinates = findCoordinates(vh, index)

            nucl = Nucleotide()
            nucl.direction = direction
            nucl.index = index
            nucl.position = [coordinates[1], coordinates[1 + direction]] #side-chain, bckbone position
            nucl.strand = strand
            nucl.vh = vh
            nucleotides_list.append(nucl)
        strand_list.append(nucleotides_list)
    return(strand_list)

def listOflistsOfNucleotides(oligos_array):
    nucl_list_list = []
    for oligo in oligos_array:
        strand_list = populateNucleotides(oligo)
        nucl_list_list.append(strand_list)
    return(nucl_list_list)

def generateVectorsandQuaternions(oligos_array):
    '''
    Given an array of oligos, fills in nucleotide attributes
    '''
    nucl_list_list = listOflistsOfNucleotides(oligos_array)

    for o, oligo in enumerate(nucl_list_list):
        for s, strand in enumerate(oligo):
            [axis_0, backbone_0] = nucl_list_list[0][0][0].position
            [axis_1, backbone_1] = nucl_list_list[0][0][1].position

            base_vector_0 = axis_0 - backbone_0
            backbone_vector_0 = backbone_1 - backbone_0
            aux_vector_a_0 = np.cross(base_vector_0, backbone_vector_0)
            aux_vector_b_0 = np.cross(aux_vector_a_0, base_vector_0)
            # return 3 orthogonal vectors in nucleotide, for quaternion
            vect_list_0 = (base_vector_0/np.linalg.norm(base_vector_0), \
                         aux_vector_a_0/np.linalg.norm(aux_vector_a_0), \
                         aux_vector_b_0/np.linalg.norm(aux_vector_b_0))

            for n, nucl in enumerate(nucl_list_list[o][s]):
                [axis_1, backbone_1] = nucl_list_list[o][s][n].position
                if n < len(nucl_list_list[o][s]) - 1:
                    [axis_2, backbone_2] = nucl_list_list[o][s][n + 1].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = backbone_2 - backbone_1
                elif n == len(nucl_list_list[o][s]) - 1:
                    [axis_2, backbone_2] = nucl_list_list[o][s][n - 1].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = - (backbone_2 - backbone_1)

                aux_vector_a_1 = np.cross(base_vector_1, backbone_vector_1)
                aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
                vect_list_1 = (base_vector_1+np.array([0.00001,0,0])/np.linalg.norm(base_vector_1+np.array([0.00001,0,0])), \
                             aux_vector_a_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_a_1+np.array([0.00001,0,0])), \
                             aux_vector_b_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_b_1+np.array([0.00001,0,0])))

                nucl.vectors = vect_list_1
                nucl.global_pts = [backbone_1, axis_1, aux_vector_a_1 + backbone_1]
                nucl_quaternion = vTools.systemQuaternion(vect_list_0, vect_list_1)
                nucl.quaternion = [nucl_quaternion.w, \
                                  nucl_quaternion.x, \
                                  nucl_quaternion.y, \
                                  nucl_quaternion.z]

    return(nucl_list_list)

# functions needed for relaxation
def distanceBetweenVhs(vh1, index1, vh2, index2):
    '''
    Given 2 points(vh, index), calculates the
    Euclian distance between them
    '''
    pos1 = findCoordinates(vh1, index1)[1]
    pos2 = findCoordinates(vh2, index2)[1]
    distance = np.linalg.norm(pos1 - pos2)
    return(distance)

def connection3p(strand):
    '''
    Given a strand, returns the vhelix to which the 3p end
    connects to, if the distance is not too far
    '''
    if strand.connection3p() != None:
            vh1 = strand.idNum()
            index1 = strand.connection3p().idx5Prime()
            vh2 = strand.connection3p().idNum()
            index2 = strand.connection3p().connection5p().idx3Prime()
            distance = distanceBetweenVhs(vh1, index1, vh2, index2)
            if distance < 10.0:
                return(vh2)

def connection5p(strand):
    '''
    Given a strand, returns the vhelix to which the 5p end
    connects to, if the distance is not too far
    '''
    if strand.connection5p() != None:
            vh1 = strand.idNum()
            index1 = strand.connection5p().idx3Prime()
            vh2 = strand.connection5p().idNum()
            index2 = strand.connection5p().connection3p().idx5Prime()
            distance = distanceBetweenVhs(vh1, index1, vh2, index2)
            if distance < 10.0:
                return(vh2)

def calculateConnections(vh):
    '''
    Given a vh number, returns the set of neighboring vhs,
    where a neighbor has a staple connection with vh
    and is closer than dist = 10.0
    '''
    staple_strandSet = part.getStrandSets(vh)[not(vh % 2)]
    connections = set()
    for strand in staple_strandSet:
        if connection3p(strand) != None:
            connections.add(connection3p(strand))
        if connection5p(strand) != None:
            connections.add(connection5p(strand))
    return(connections)

def separateOrigamiParts(part):
    '''
    Separates the origami 'part' into bodies
    by evaluating if a vh is connected to others
    see: 'calculateConnections'
    '''
    vhs_list = list(part.getIdNums())
    bodies = []

    for vh in vhs_list: #loop over the vhs of part
        body_index = None
        vh_connections = calculateConnections(vh)
        for b, body in enumerate(bodies): #loop over potential bodies already seen
            if vh in body: #this vh is part of a known body
                body_index = b
                break
            elif vh_connections.intersection(body) != set():
            #one of the connections belong to an existing body
                body_index = b
                break
        if body_index == None: # not vh nor its connections are in known bodies
            body_index = len(bodies)
            bodies.append(set())

        bodies[body_index].add(vh)
        bodies[body_index].update(vh_connections)
    return(bodies)

def populateBodiesNuclAndVhs(nucleotides_list_of_list):
    '''
    Given a list of oligos, each composed of a list of strands
    each composed of a list of nucleotides, find out which body
    each nucleotide belongs to and add it to the corresponding
    body set.
    '''
    vhs_of_body = separateOrigamiParts(part)
    vhs_of_body_list = list(vhs_of_body)
    num_bodies = len(vhs_of_body_list)
    bodies = [Body() for i in range(num_bodies)]

    for chain in nucleotides_list_of_list:
        for strand in chain:
            for nucl in strand:
                # seach in each body for this nucleotides' vh, assign to body
                body_id = [i for i in range(num_bodies) if nucl.vh in vhs_of_body[i]][0]
                nucl.body = body_id
                bodies[body_id].add_nucleotide(nucl)
                bodies[body_id].add_vh(nucl.vh)
    return(bodies)

def populateBody(nucleotides_list_of_list):
    '''
    Given a list of oligos, each composed of a list of strands
    each composed of a list of nucleotides,
    first populate each body's nucleotide and Vh and
    then calculate the other attributes.
    '''

    bodies = populateBodiesNuclAndVhs(nucleotides_list_of_list)

    for body in bodies:
        positions = [nucl.position[1] for nucl in body.nucleotides]
        body.com_position = vTools.calculateCoM(positions)
        body.moment_inertia = vTools.calculateMomentInertia(positions)
        body.com_quaternion = [1., 0., 0., 0.]
    return(bodies)

###################################
# start cadnano and read input file
###################################
# Read nucleotides data from cadnano
# Read nucleotides data from cadnano
app = cadnano.app()
doc = app.document = Document()
doc.readFile(INPUT_FILENAME);

part = doc.activePart()
oligos = part.oligos()

oligos_sorted_by_length = sorted(oligos, key=lambda x: x.length(), reverse=True)
longest_oligo = oligos_sorted_by_length[0]
staple_oligos = oligos_sorted_by_length[1:]
oligos_array = [longest_oligo] + [staple for staple in staple_oligos]
# calculate list of oligos, each as a list of nucleotides in 5'-3' direction
list_of_list_of_nucleotides = generateVectorsandQuaternions(oligos_array)

##################################
#start HOOMD code
##################################
from hoomd import *
from hoomd import md

# Start HOOMD
context.initialize("");



if RELAX:
    bodies = populateBody(list_of_list_of_nucleotides)








#
#
