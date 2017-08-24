from os import path
import numpy as np
import cadnano
from cadnano.document import Document
import vector_tools as vTools #needed for vector / quaternion calculations

####### USER DEFINED AND GLOBAL VARIABLES ######
INPUT_FILENAME = 'input/tripod.json'
RELAX = True

global_nucl_matrix = []
global_connections_pairs = []
###############################################

##################################
#  Helper functions for cadnano  #
##################################
class Nucleotide:
    '''
    Fixed attributes of a nucleotide
    Attributes: index, position, strand, vh
    '''
    def __init__(self):
        self.direction  = 1    #1 is fwd, -1 is reverse
        self.global_pts = None # backbone, sidechain and aux points in global frame
        self.index      = None # z position in cadnano's unit
        self.body       = None # body this nucl belongs to (for relaxation)
        self.ext_conn   = None # vh and id of nucleotide connected in other body
        self.position   = None # positions of axis (sidechain) and backbone ptcls
        self.quaternion = None # quaternion orientation for this nucleotide
        self.strand	    = None # strand #
        self.vectors    = None # orth vectors in the body ref. frame for quat calc
        self.vh  	    = None # virtual helix this nucleotide belongs to

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
    (strand, strand direction, virtual_helix_id, index)
    for each nucleotide. List is oriented 5' to 3'
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

def oligosArray(active_part):
    '''
    Given a origami part
    Return an array with all oligos in part sorted by length
    '''
    oligos = active_part.oligos()
    oligos_sorted_by_length = sorted(oligos, key=lambda x: x.length(), reverse=True)
    return(oligos_sorted_by_length)

def findCoordinates(vh, index):
    '''
    Given a vh and a index, returns (x,y,z)
    for the sidechain pts and backbones fwd and rev
    '''
    axis_pts = part.getCoordinates(vh)[0][index]
    fwd_pts = part.getCoordinates(vh)[1][index]
    rev_pts = part.getCoordinates(vh)[2][index]
    return [fwd_pts, axis_pts, rev_pts]

def populateGlobalNuclMatrix(active_part):
    '''
    Creates an empty matrix of len = vh_length x index_length
    to be populated with all nucleotides in part
    This will be the global reference to any nucleotide
    via global_nucl_matrix[vh][index][rev or fwd]
    '''
    global global_nucl_matrix
    vhs_length = len(list(active_part.getIdNums()))
    bases_length = active_part.getVirtualHelix(0).getSize()
    global_nucl_matrix = [[[[] for k in range(2)] for i in range(bases_length)] for j in range(vhs_length)]

def populateBasicNucleotideAttributes(oligo, active_part):
    '''
    Given an oligo, returns a list of strands,
    each containing the pointers ([vh][index][is_fwd]) to the
    nucleotides making up such strand and populate nucleotides matrix
    with basic attributes (direction, index, position, strand, vh)
    '''
    global global_nucl_matrix
    if global_nucl_matrix == []:
        populateGlobalNuclMatrix(active_part) #start by pre-populating global var

    strand_list = []
    for helper_strands in oligoHelperList(oligo):
        nucleotides_list = []
        for [strand, direction, vh, index] in helper_strands:
            coordinates = findCoordinates(vh, index)
            is_fwd = int((direction + 1)/2)
            nucl = Nucleotide()
            nucl.direction = direction
            nucl.index = index
            nucl.position = [coordinates[1], coordinates[1 + direction]] #side-chain, bckbone position
            nucl.strand = strand
            nucl.vh = vh
            global_nucl_matrix[vh][index][is_fwd] = nucl

            nucleotides_list.append([vh, index, is_fwd])
        strand_list.append(nucleotides_list)
    return(strand_list)

def createOligosList(oligos_array, active_part):
    '''
    Given an array of oligos in part,
    returns a list of oligos, each containing
    a list of strands, each containing a list of nucleotides
    making up the part.
    '''
    global global_nucl_matrix
    oligosList = []
    for oligo in oligos_array:
        strand_list = populateBasicNucleotideAttributes(oligo, active_part)
        oligosList.append(strand_list)
    return(oligosList)

def translateIndices(oligos_helper_list, i, j, k):
    '''
    Oligos helper list gives a tuple list of pointers
    of the type [vh, index, is_fwd] which are needed to
    reference nucleotides in global_nucl_matrix
    this function translates oligo_helper_list into
    indices to be used by global_nucl_matrix
    '''
    [vh, index, is_fwd] = oligos_helper_list[i][j][k]
    return([vh, index, is_fwd])

def populateAllNucleotideAttributes(oligos_helper_list, active_part):
    '''
    Given an array of oligos, fills the remaining nucleotide attributes:
    vectors, quaternion, global_pts
    '''
    global global_nucl_matrix

    for o, oligo in enumerate(oligos_helper_list):
        for s, strand in enumerate(oligo):

            [vh_0, index_0, is_fwd_0] = translateIndices(oligos_helper_list, 0, 0, 0)
            [vh_1, index_1, is_fwd_1] = translateIndices(oligos_helper_list, 0, 0, 1)

            [axis_0, backbone_0] = global_nucl_matrix[vh_0][index_0][is_fwd_0].position
            [axis_1, backbone_1] = global_nucl_matrix[vh_1][index_1][is_fwd_1].position

            base_vector_0 = axis_0 - backbone_0
            backbone_vector_0 = backbone_1 - backbone_0
            aux_vector_a_0 = np.cross(base_vector_0, backbone_vector_0)
            aux_vector_b_0 = np.cross(aux_vector_a_0, base_vector_0)
            # return 3 orthogonal vectors in nucleotide, for quaternion
            vect_list_0 = (base_vector_0/np.linalg.norm(base_vector_0), \
                         aux_vector_a_0/np.linalg.norm(aux_vector_a_0), \
                         aux_vector_b_0/np.linalg.norm(aux_vector_b_0))

            for n, nucl in enumerate(oligos_helper_list[o][s]):
                [vh_1, index_1, is_fwd_1] = translateIndices(oligos_helper_list, o, s, n)
                [axis_1, backbone_1] = global_nucl_matrix[vh_1][index_1][is_fwd_1].position

                if n < len(oligos_helper_list[o][s]) - 1:
                    [vh_2, index_2, is_fwd_2] = translateIndices(oligos_helper_list, o, s, n + 1)
                    [axis_2, backbone_2] = global_nucl_matrix[vh_2][index_2][is_fwd_2].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = backbone_2 - backbone_1
                elif n == len(oligos_helper_list[o][s]) - 1:
                    [vh_2, index_2, is_fwd_2] = translateIndices(oligos_helper_list, o, s, n - 1)
                    [axis_2, backbone_2] = global_nucl_matrix[vh_2][index_2][is_fwd_2].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = - (backbone_2 - backbone_1)

                aux_vector_a_1 = np.cross(base_vector_1, backbone_vector_1)
                aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
                vect_list_1 = (base_vector_1+np.array([0.00001,0,0])/np.linalg.norm(base_vector_1+np.array([0.00001,0,0])), \
                             aux_vector_a_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_a_1+np.array([0.00001,0,0])), \
                             aux_vector_b_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_b_1+np.array([0.00001,0,0])))

                this_nucl = global_nucl_matrix[vh_1][index_1][is_fwd_1]
                this_nucl.vectors = vect_list_1
                this_nucl.global_pts = [backbone_1, axis_1, aux_vector_a_1 + backbone_1]
                nucl_quaternion = vTools.systemQuaternion(vect_list_0, vect_list_1)
                this_nucl.quaternion = [nucl_quaternion.w, \
                                        nucl_quaternion.x, \
                                        nucl_quaternion.y, \
                                        nucl_quaternion.z]
    return()

###################################
# Functions needed for relaxation #
###################################
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
    global global_connections_pairs
    if strand.connection3p() != None:
            vh1 = strand.idNum()
            index2 = strand.connection3p().idx5Prime()
            vh2 = strand.connection3p().idNum()
            index1 = strand.connection3p().connection5p().idx3Prime()
            distance = distanceBetweenVhs(vh1, index1, vh2, index2)
            if distance < 6.0:
                return(vh2)
            else:
                global_connections_pairs.append([[vh1, index1], [vh2, index2]])

def connection5p(strand):
    '''
    Given a strand, returns the vhelix to which the 5p end
    connects to, if the distance is not too far
    '''
    global global_connections_pairs
    if strand.connection5p() != None:
            vh1 = strand.idNum()
            index2 = strand.connection5p().idx3Prime()
            vh2 = strand.connection5p().idNum()
            index1 = strand.connection5p().connection3p().idx5Prime()
            distance = distanceBetweenVhs(vh1, index1, vh2, index2)
            if distance < 6.0:
                return(vh2)
            else:
                global_connections_pairs.append([[vh1, index1], [vh2, index2]])

def calculateConnections(vh):
    '''
    Given a vh number, returns the set of neighboring vhs,
    where a neighbor has a *staple* connection with vh
    that is closer than dist = 10.0
    '''
    staple_strandSet = part.getStrandSets(vh)[not(vh % 2)] #staple strands only
    scaffold_strandSet = part.getStrandSets(vh)[(vh % 2)]
    connections = set()
    for strand in staple_strandSet:
        if connection3p(strand) != None:
            connections.add(connection3p(strand))
        if connection5p(strand) != None:
            connections.add(connection5p(strand))
    for strand in scaffold_strandSet:
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
            #one of the connections belongs to an known body
                body_index = b
                break
        if body_index == None: # not vh nor its connections are in known bodies
            body_index = len(bodies)
            bodies.append(set())

        bodies[body_index].add(vh)
        bodies[body_index].update(vh_connections)
    return(bodies)

def populateBodiesNuclAndVhs(oligos_helper_list):
    '''
    Given a list of oligos, each composed of a list of strands
    each composed of a list of nucleotides, find out which body
    each nucleotide belongs to and add it to the corresponding
    body set.
    '''
    global global_nucl_matrix

    vhs_of_body = separateOrigamiParts(part)
    vhs_of_body_list = list(vhs_of_body)
    num_bodies = len(vhs_of_body_list)
    bodies = [Body() for i in range(num_bodies)]

    for o, oligo in enumerate(oligos_helper_list):
        for s, strand in enumerate(oligo):
            for n, nucl in enumerate(strand):
                # seach in each body for this nucleotides' vh, assign to body
                [vh_1, index_1, is_fwd_1] = translateIndices(oligos_helper_list, o, s, n)
                this_nucl = global_nucl_matrix[vh_1][index_1][is_fwd_1]
                body_id = [i for i in range(num_bodies) if this_nucl.vh in vhs_of_body[i]][0]
                this_nucl.body = body_id
                bodies[body_id].add_nucleotide(this_nucl)
                bodies[body_id].add_vh(this_nucl.vh)
    return(bodies)

def populateBody(oligos_helper_list):
    '''
    Given a list of oligos, each composed of a list of strands
    each composed of a list of nucleotides,
    first populate each body's nucleotide and Vh and
    then calculate the other attributes.
    '''
    bodies = populateBodiesNuclAndVhs(oligos_helper_list)

    for body in bodies:
        positions = [nucl.position[1] for nucl in body.nucleotides]
        body.com_position = vTools.calculateCoM(positions)
        body.moment_inertia = vTools.calculateMomentInertia(positions)
        body.com_quaternion = [1., 0., 0., 0.]
    return(bodies)


###################################
# start cadnano and read input file
###################################

app = cadnano.app()
doc = app.document = Document()
doc.readFile(INPUT_FILENAME);
part = doc.activePart()
oligos_array = oligosArray(part)

oligos_helper_list = createOligosList(oligos_array, part)
populateAllNucleotideAttributes(oligos_helper_list, part)

bodies = populateBody(oligos_helper_list)

bodies[5].vhs


##################################
#start HOOMD code
##################################
# from hoomd import *
# from hoomd import md
#
# # Start HOOMD
# context.initialize("");
#
# num_rigid_bodies = len(bodies)
#
# particle_positions = [nucl.positions for ]














#
#
