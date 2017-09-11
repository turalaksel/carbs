import numpy as np
import cadnano
import vectortools
import functools

from cadnano.document import Document
from hoomd import *
from hoomd import md

class Origami:
    '''
    Parent DNA origami model class
    '''

    def __init__(self):
        #Body variables
        self.rigid_bodies                  = []
        self.soft_bodies                   = []


        self.rigid_body_nucleotide_types   = []
        self.rigid_body_nucleotides        = []

        self.soft_body_nucleotide_types    = []
        self.soft_body_nucleotides         = []

        self.num_rigid_bodies              = 0
        self.num_soft_bodies               = 0

        #Soft connections
        self.inter_rigid_body_connections  = []
        self.inter_nucleotide__connections = []

        #Cadnano parameters
        self.part                   = None
        self.oligos                 = None
        self.strands                = None
        self.num_vhs                = None
        self.nucleotide_list        = None
        self.nucleotide_matrix      = None

        self.nucleotide_type_list   = None
        self.nucleotide_type_matrix = None

        self.crossovers             = None
        self.vh_vh_crossovers       = None
        self.long_range_connections = {}
        self.short_range_connections= {}
        self.soft_connections       = {}

        #Distance constraints
        self.crossover_distance      = 2.0   #Distance in Angstrom

    def parse_soft_connections(self):
        self.inter_rigid_body_connections = set()
        self.inter_nucleotide_connections = set()

        for pointer_1, pointer_2 in self.soft_connections.items():
            vh_1,index_1,is_fwd_1 = pointer_1
            vh_2,index_2,is_fwd_2 = pointer_2

            nucleotide_1 = self.nucleotide_matrix[vh_1][index_1][is_fwd_1]
            nucleotide_2 = self.nucleotide_matrix[vh_2][index_2][is_fwd_2]

            #Get body numbers
            body_num_1   = nucleotide_1.body_num
            body_num_2   = nucleotide_2.body_num

            #Skip the skip nucleotides
            if nucleotide_1.skip or nucleotide_2.skip:
                continue

            #If both nucleotides are part of the rigid bodies, add the connection to rigid body connections
            if nucleotide_1.body.type and nucleotide_2.body.type and body_num_1 != body_num_2 and not (body_num_2,body_num_1) in self.inter_rigid_body_connections:
                self.inter_rigid_body_connections.add((body_num_1,body_num_2))


            #Add nucleotide-nucleotide connections to list
            nucleotide_sim_num_1 = nucleotide_1.simulation_nucleotide_num
            nucleotide_sim_num_2 = nucleotide_2.simulation_nucleotide_num

            if not (nucleotide_sim_num_2,nucleotide_sim_num_1) in self.inter_nucleotide_connections:
                self.inter_nucleotide_connections.add((nucleotide_sim_num_1,nucleotide_sim_num_2))

    def incorporate_skips(self):
        '''
        Incorporate skips in the model
        '''
        #1. Determine the skip boundaries
        self.skip_boundaries  = []
        for vh in range(self.num_vhs):
            for is_fwd in [0,1]:

                skip_begin_found = False
                skip_begin = None
                skip_end   = None
                for idx in range(len(self.skip_matrix[vh])):
                    if self.skip_matrix[vh][idx][is_fwd] and not skip_begin_found:
                        skip_begin_found = True
                        skip_begin = idx
                    elif not self.skip_matrix[vh][idx][is_fwd] and skip_begin_found:
                        skip_begin_found = False
                        skip_end = idx - 1
                        self.skip_boundaries.append([vh,is_fwd,skip_begin,skip_end])

        #2. Connect the beginning and end of
        for vh,is_fwd,idx_begin,idx_end in self.skip_boundaries:
            pointer_1 = (vh,idx_begin-1,1)
            pointer_2 = (vh,idx_end+1,1)

            if not is_fwd:
                pointer_1 = (vh,idx_end+1,0)
                pointer_2 = (vh,idx_begin-1,0)

            self.short_range_connections[pointer_1] = pointer_2

    def assign_nucleotide_types(self):
        '''
        Build the nucleotide network for an origami design
        '''
        self.nucleotide_type_list = []
        for vh in range(len(self.nucleotide_matrix)):
            for idx in range(len(self.nucleotide_matrix[vh])):

                #Determine the nucleotide type ssNucleotide vs dsNucletide and nucleotide connections
                current_nucleotide_rev = self.nucleotide_matrix[vh][idx][0]
                current_nucleotide_fwd = self.nucleotide_matrix[vh][idx][1]

                if not current_nucleotide_fwd == None and not current_nucleotide_rev == None:
                    ds_nucleotide = DSNucleotide()
                    ds_nucleotide.fwd_nucleotide   = current_nucleotide_fwd
                    ds_nucleotide.rev_nucleotide   = current_nucleotide_rev
                    ds_nucleotide.type             = 1
                    ds_nucleotide.skip             = current_nucleotide_fwd.skip
                    self.nucleotide_type_matrix[vh][idx] = ds_nucleotide
                    self.nucleotide_type_list.append(ds_nucleotide)

                elif not current_nucleotide_fwd == None:
                    ss_nucleotide                  = SSNucleotide()
                    ss_nucleotide.nucleotide       = current_nucleotide_fwd
                    ss_nucleotide.type             = 0
                    ss_nucleotide.skip             = current_nucleotide_fwd.skip

                    self.nucleotide_type_matrix[vh][idx] = ss_nucleotide
                    self.nucleotide_type_list.append(ss_nucleotide)

                elif not current_nucleotide_rev == None:

                    ss_nucleotide                  = SSNucleotide()
                    ss_nucleotide.nucleotide       = current_nucleotide_rev
                    ss_nucleotide.type             = 0
                    ss_nucleotide.skip             = current_nucleotide_rev.skip

                    self.nucleotide_type_matrix[vh][idx] = ss_nucleotide
                    self.nucleotide_type_list.append(ss_nucleotide)

    def assign_nucleotide_connections(self):
        '''
        Assign nucleotide connections
        '''
        #1. Add the base-stacking interactions
        for vh in range(len(self.nucleotide_type_matrix)):
            for idx in range(len(self.nucleotide_type_matrix[vh])):

                if self.nucleotide_type_matrix[vh][idx] == None or self.nucleotide_type_matrix[vh][idx].skip:
                    continue

                self.nucleotide_type_matrix[vh][idx].rigid_connections = []
                self.nucleotide_type_matrix[vh][idx].soft_connections  = []

                #Get the type for nucleotide
                type_1 = self.nucleotide_type_matrix[vh][idx].type

                #Pointer 1
                pointer1_fwd = (vh,idx,1)
                pointer1_rev = (vh,idx,0)

                if  idx+1 < len(self.nucleotide_type_matrix[vh]) and not self.nucleotide_type_matrix[vh][idx+1] == None and not self.nucleotide_type_matrix[vh][idx+1].skip:

                    type_2 = self.nucleotide_type_matrix[vh][idx+1].type

                    #If both types are dsNucleotide make the connection rigid
                    if type_1*type_2:
                        self.nucleotide_type_matrix[vh][idx].rigid_connections.append(self.nucleotide_type_matrix[vh][idx+1])
                    else:
                        self.nucleotide_type_matrix[vh][idx].soft_connections.append(self.nucleotide_type_matrix[vh][idx+1])

                        pointer2_fwd = (vh,idx+1,1)
                        pointer2_rev = (vh,idx+1,0)

                        if  pointer1_fwd in  self.short_range_connections.keys() and self.short_range_connections[pointer1_fwd] == pointer2_fwd:
                            self.soft_connections[pointer1_fwd] = pointer2_fwd

                        elif self.short_range_connections[pointer2_rev] == pointer1_rev:
                            self.soft_connections[pointer2_rev] = pointer1_rev

                if idx-1 >= 0 and not self.nucleotide_type_matrix[vh][idx-1] == None and not self.nucleotide_type_matrix[vh][idx-1].skip:
                    type_2 = self.nucleotide_type_matrix[vh][idx-1].type
                    if type_1*type_2:
                        self.nucleotide_type_matrix[vh][idx].rigid_connections.append(self.nucleotide_type_matrix[vh][idx-1])
                    else:
                        self.nucleotide_type_matrix[vh][idx].soft_connections.append(self.nucleotide_type_matrix[vh][idx-1])

                        pointer2_fwd = (vh,idx-1,1)
                        pointer2_rev = (vh,idx-1,0)

                        if pointer2_fwd in  self.short_range_connections.keys() and self.short_range_connections[pointer2_fwd] == pointer1_fwd:
                            self.soft_connections[pointer2_fwd] = pointer1_fwd

                        elif self.short_range_connections[pointer1_rev] == pointer2_rev:
                            self.soft_connections[pointer1_rev] = pointer2_rev

        #2. Add short range connections that are not adjacent on sequence due to skips
        for pointer1,pointer2 in self.short_range_connections.items():
            vh1,idx1,is_fwd1 = pointer1
            vh2,idx2,is_fwd2 = pointer2

            #If the bases are not adjacent in sequence, add the connections to soft connections
            if abs(idx1-idx2) > 1:
                #Add the connections first in nucleotide type matrix
                self.nucleotide_type_matrix[vh1][idx1].soft_connections.append(self.nucleotide_type_matrix[vh2][idx2])
                self.nucleotide_type_matrix[vh1][idx1].soft_connections.append(self.nucleotide_type_matrix[vh2][idx2])
                self.soft_connections[pointer1] = pointer2
                print(pointer1,pointer2)

        #3. Add the crossover connections
        for pointer_1, pointer_2 in self.crossovers.items():
            (vh_1, index_1, is_fwd_1) = pointer_1
            (vh_2, index_2, is_fwd_2) = pointer_2

            type_1 = self.nucleotide_type_matrix[vh_1][index_1].type
            type_2 = self.nucleotide_type_matrix[vh_2][index_2].type

            if self.vh_vh_crossovers[vh_1][vh_2] > 1 and type_1*type_2:
                self.nucleotide_type_matrix[vh_1][index_1].rigid_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
                self.nucleotide_type_matrix[vh_2][index_2].rigid_connections.append(self.nucleotide_type_matrix[vh_1][index_1])
            else:
                self.nucleotide_type_matrix[vh_1][index_1].soft_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
                self.nucleotide_type_matrix[vh_2][index_2].soft_connections.append(self.nucleotide_type_matrix[vh_1][index_1])

                #Add the connection to soft connection list
                self.soft_connections[pointer_1] = pointer_2

        #4. Add long range connections
        for pointer_1, pointer_2 in self.long_range_connections.items():
            (vh_1, index_1, is_fwd_1) = pointer_1
            (vh_2, index_2, is_fwd_2) = pointer_2

            self.nucleotide_type_matrix[vh_1][index_1].soft_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
            self.nucleotide_type_matrix[vh_2][index_2].soft_connections.append(self.nucleotide_type_matrix[vh_1][index_1])

            #Add the connection to soft connection list
            self.soft_connections[pointer_1] = pointer_2

    def get_connections(self):
        '''
        Given a vh number, returns the set of neighboring vhs,
        where a neighbor has a *staple* connection with vh
        that is closer than dist = 10.0
        '''
        self.crossovers             = {}
        self.long_range_connections = {}
        for vh in range(self.num_vhs):
            staple_strandSet   = self.part.getStrandSets(vh)[not(vh % 2)] #staple strands only
            scaffold_strandSet = self.part.getStrandSets(vh)[(vh % 2)]

            for strand in staple_strandSet:
                self.connection3p(strand)
            for strand in scaffold_strandSet:
                self.connection3p(strand)

    def dfs(self, start_nucleotide_type):
        '''
        Depth-first-search graph traverse algorithm to find connected components
        '''
        visited, stack = set(), [start_nucleotide_type]
        while stack:
            new_nucleotide_type = stack.pop()
            if new_nucleotide_type not in visited:
                visited.add(new_nucleotide_type)
                new_nucleotide_type.visited = True
                stack.extend(set(new_nucleotide_type.rigid_connections) - visited)
        return visited

    def cluster_into_bodies(self):
        '''
        Cluster the DNA origami structure into body clusters
        '''

        self.num_bodies            = 0
        self.body_nucleotide_types = []
        self.body_nucleotides      = []
        self.rigid_bodies          = []
        self.soft_bodies           = []

        #1. Identify the clusters using dfs
        for nucleotide_type in self.nucleotide_type_list:

            if nucleotide_type.visited or nucleotide_type.skip:
                continue

            #Get the nucleotide types for each cluster
            self.body_nucleotide_types.append([])
            self.body_nucleotide_types[self.num_bodies] = list(self.dfs(nucleotide_type))

            #Check if the cluster is a rigid body
            rigid_body = functools.reduce(lambda x,y:x*y,[nucleotide_type.type for nucleotide_type in self.body_nucleotide_types[self.num_bodies]])

            #Create a new Body object
            new_body = Body()
            new_body.nucleotide_types = self.body_nucleotide_types[self.num_bodies]
            new_body.type             = rigid_body

            #If the body is rigid add to rigid body collection
            if rigid_body:
                self.rigid_bodies.append(new_body)
            else:
                self.soft_bodies.append(new_body)

            self.num_bodies += 1

        #2. Update soft body nucleotide body position numbers and body numbers
        nucleotide_number = 0
        for i in range(len(self.soft_bodies)):
            soft_body               = self.soft_bodies[i]
            soft_body.nucleotides   = []
            nucleotide_type         = soft_body.nucleotide_types[0]

            #Update the nucleotide number
            nucleotide_type.nucleotide.simulation_nucleotide_num = nucleotide_number

            #Assign soft body to nucleotide
            nucleotide_type.nucleotide.body = soft_body

            soft_body.nucleotides  += [nucleotide_type.nucleotide]
            self.body_nucleotides  += [nucleotide_type.nucleotide]

            nucleotide_number      += 1

            #Initialize soft bodies
            soft_body.initialize()

        #3. Update rigid body nucleotide body position numbers and body numbers
        for i in range(len(self.rigid_bodies)):
            rigid_body             = self.rigid_bodies[i]
            rigid_body.nucleotides = []
            for nucleotide_type in rigid_body.nucleotide_types:
                fwd_nucleotide = nucleotide_type.fwd_nucleotide
                rev_nucleotide = nucleotide_type.rev_nucleotide

                fwd_nucleotide.simulation_nucleotide_num = nucleotide_number
                rev_nucleotide.simulation_nucleotide_num = nucleotide_number+1

                #Assign the rigid bodies to nucleotides
                fwd_nucleotide.body = rigid_body
                rev_nucleotide.body = rigid_body

                nucleotide_number    += 2

                rigid_body.nucleotides+= [fwd_nucleotide,rev_nucleotide]
                self.body_nucleotides += [fwd_nucleotide,rev_nucleotide]

            #Initialize rigid bodies
            rigid_body.initialize()

    def parse_oligo(self, oligo):
        '''
        Given an oligo, returns the following list of useful properties
        (strand, strand direction, virtual_helix_id, index)
        for each nucleotide. List is oriented 5' to 3'
        '''

        generator         = oligo.strand5p().generator3pStrand()
        oligo_helper_list = []
        for strand in generator:
            strand_helper_list = []
            vh        = strand.idNum()
            index_5p  = strand.idx5Prime()
            index_3p  = strand.idx3Prime()
            direction = (-1 + 2*strand.isForward()) #-1 if backwards, 1 if fwd
            for i in range(index_5p, index_3p + direction, direction):
                strand_helper_list.append([strand, direction, vh, i])
            oligo_helper_list.append(strand_helper_list)
        return oligo_helper_list

    def find_skips(self):
        '''
        Identify all the skips in the structure
        '''
        for oligo in self.oligos:
            generator  = oligo.strand5p().generator3pStrand()
            for strand in generator:
                vh = strand.idNum()
                index_5p  = strand.idx5Prime()
                index_3p  = strand.idx3Prime()
                direction = (-1 + 2*strand.isForward()) #-1 if backwards, 1 if fwd
                for idx in range(index_5p, index_3p + direction, direction):
                    if strand.hasInsertionAt(idx) and strand.insertionLengthBetweenIdxs(idx,idx) == -1:
                        self.skip_matrix[vh][idx][int(strand.isForward())] = True

    def list_oligos(self):
        '''
        Given a origami part
        Return an array with all oligos in part sorted by length
        '''
        self.oligos = self.part.oligos()
        self.oligos = sorted(self.oligos, key=lambda x: x.length(), reverse=True)

        return self.oligos

    def get_coordinates(self, vh, index):
        '''
        Given a vh and a index, returns (x,y,z)
        for the sidechain pts and backbones fwd and rev
        '''
        axis_pts = self.part.getCoordinates(vh)[0][index]
        fwd_pts  = self.part.getCoordinates(vh)[1][index]
        rev_pts  = self.part.getCoordinates(vh)[2][index]

        return [rev_pts, fwd_pts, axis_pts]

    def initialize_nucleotide_matrix(self):
        '''
        nucleotide_matrix is a 2D list of nucleotides at (vh, idx).
        nucleotide_type_matrix is whether a nucleotide is ssDNA or dsDNA
        vh_vh_crossovers is whether a vh_1 crosses over to vh_2
        skip_matrix is 2D list of possible skips at nucleotide (vh, idx)
        '''
        self.num_vhs = len(list(self.part.getIdNums()))
        num_bases    = self.part.getVirtualHelix(0).getSize()

        self.nucleotide_matrix       = [[[None,None]  for idx in range(num_bases)] for vh in range(self.num_vhs)]
        self.nucleotide_type_matrix  = [[ None  for idx in range(num_bases)] for vh in range(self.num_vhs)]
        self.vh_vh_crossovers        = [[0  for vh in range(self.num_vhs)] for vh in range(self.num_vhs)]
        self.skip_matrix             = [[[False,False]  for idx in range(num_bases)] for vh in range(self.num_vhs)]


    def connection5p(self,strand):
        '''
        Given a strand, find its 5' crossovers. If connection distance is within
        self.crossover_distance cutoff, add to vh_vh_crossovers list.
        Else, add to 'long_range_connections' list
        '''
        if strand.connection5p() != None:
            vh_1     = strand.idNum()
            index_1  = strand.idx5Prime()
            is_fwd_1 = int(strand.isForward())

            vh_2     = strand.connection5p().idNum()
            index_2  = strand.connection5p().idx3Prime()
            is_fwd_2 = int(strand.connection5p().isForward())

            conn_pointer_1 = (vh_1, index_1, is_fwd_1)
            conn_pointer_2 = (vh_2, index_2, is_fwd_2)

            distance = self.distance_between_vhs(vh_1, index_1, is_fwd_1, vh_2, index_2, is_fwd_2)

            if distance < self.crossover_distance:
                self.crossovers[conn_pointer_1]    =  conn_pointer_2
                self.vh_vh_crossovers[vh_1][vh_2] += 1
                self.vh_vh_crossovers[vh_2][vh_1] += 1
            else:
                self.long_range_connections[conn_pointer_1] = conn_pointer_2


    def connection3p(self,strand):
        '''
        Given a strand, find its 3' crossovers. If connection distance is within
        self.crossover_distance cutoff, add to vh_vh_crossovers list.
        Else, add to 'long_range_connections' list
        '''
        if strand.connection3p() != None:
            vh_1 = strand.idNum()
            index_1 = strand.idx3Prime()
            is_fwd_1 = int(strand.isForward())

            vh_2 = strand.connection3p().idNum()
            index_2 = strand.connection3p().idx5Prime()
            is_fwd_2 = int(strand.connection3p().isForward())

            conn_pointer_1 = (vh_1, index_1, is_fwd_1)
            conn_pointer_2 = (vh_2, index_2, is_fwd_2)

            distance = self.distance_between_vhs(vh_1, index_1, is_fwd_1, vh_2, index_2, is_fwd_2)

            if distance <self.crossover_distance:
                self.crossovers[conn_pointer_1]    =  conn_pointer_2
                self.vh_vh_crossovers[vh_1][vh_2] += 1
                self.vh_vh_crossovers[vh_2][vh_1] += 1
            else:
                self.long_range_connections[conn_pointer_1] = conn_pointer_2

    def get_nucleotide(self,pointers):
        '''
        Given a tuple of pointers in the form [vh, index, is_fwd],
        Returns the global nucleotide referent to the pointers
        '''
        [vh, index, is_fwd] = pointers
        return self.nucleotide_matrix[vh][index][is_fwd]

    def oligo_to_strands_nucleotides(self, oligo):
        '''
        Given an oligo, returns a list of strands,
        each containing the pointers ([vh][index][is_fwd]) to the
        nucleotides making up such strand and populate nucleotides matrix
        with basic attributes (direction, index, position, strand, vh)
        '''
        if self.nucleotide_matrix == None:
            self.initialize_nucleotide_matrix()                             #start pre-populating global var

        strand_list = []
        for helper_strands in self.parse_oligo(oligo):
            nucleotides_list = []
            for i in range(len(helper_strands)):
                #Current nucleotide
                strand, direction, vh, index = helper_strands[i]

                if i+1 < len(helper_strands):
                    strand_next, direction_next,vh_next,index_next = helper_strands[i+1]
                    conn_pointer_1 = (vh     ,index     , int(direction > 0     ))
                    conn_pointer_2 = (vh_next,index_next, int(direction_next > 0))

                    #Assign short range connection
                    self.short_range_connections[conn_pointer_1] = conn_pointer_2

                #Get the coordinates
                coordinates = self.get_coordinates(vh, index)

                #Get direction
                is_fwd = int(direction > 0)

                new_nucleotide = Nucleotide()
                new_nucleotide.direction         = direction
                new_nucleotide.index             = index
                new_nucleotide.position          = [coordinates[2], coordinates[is_fwd]] #Sidechain(axis) and backbone coordinates
                new_nucleotide.strand            = strand
                new_nucleotide.vh                = vh
                new_nucleotide.is_fwd            = is_fwd
                new_nucleotide.skip              = self.skip_matrix[vh][index][is_fwd]

                #Assign the nucleotide
                self.nucleotide_matrix[vh][index][is_fwd] = new_nucleotide

                #Add nucleotide to list
                self.nucleotide_list.append(new_nucleotide)

                nucleotides_list.append([vh, index, is_fwd])
            strand_list.append(nucleotides_list)
        return strand_list

    def oligos_list_to_nucleotide_info(self, i, j, k):
        '''
        Oligos helper list gives a tuple list of pointers
        of the type [vh, index, is_fwd] which are needed to
        reference nucleotides in global_nucl_matrix
        this function translates oligo_helper_list into
        indices to be used by global_nucl_matrix
        '''

        [vh, index, is_fwd] = self.oligos_list[i][j][k]
        return [vh, index, is_fwd]

    def create_oligos_list(self):
        '''
        Given an array of oligos in part, returns a list of *oligos*,
        each containing a list of *strands(), each containing a
        list of *nucleotides() making up the part.
        '''
        self.oligos_list     = []
        self.nucleotide_list = []
        for oligo in self.oligos:
            strand_list = self.oligo_to_strands_nucleotides(oligo)
            self.oligos_list.append(strand_list)
        return self.oligos_list

    def distance_between_vhs(self, vh1, index1, is_fwd1, vh2, index2, is_fwd2):
        '''
        Given 2 points(vh, index), calculates the
        Euclian distance between them
        '''
        pos1 = self.get_coordinates(vh1, index1)[is_fwd1]
        pos2 = self.get_coordinates(vh2, index2)[is_fwd2]
        distance = np.linalg.norm(pos1 - pos2)
        return distance

    def populate_nucleotide_geometries(self):
        '''
        Given an array of oligos, fills the remaining nucleotide attributes:
        vectors, quaternion, global_pts
        '''

        for o, oligo in enumerate(self.oligos_list):
            for s, strand in enumerate(oligo):

                [vh_0, index_0, is_fwd_0] = self.oligos_list_to_nucleotide_info( 0, 0, 0)
                [vh_1, index_1, is_fwd_1] = self.oligos_list_to_nucleotide_info( 0, 0, 1)

                [axis_0, backbone_0] = self.nucleotide_matrix[vh_0][index_0][is_fwd_0].position
                [axis_1, backbone_1] = self.nucleotide_matrix[vh_1][index_1][is_fwd_1].position

                #Calculate the vectors
                base_vector_0     = axis_0     - backbone_0
                backbone_vector_0 = backbone_1 - backbone_0

                aux_vector_a_0 = np.cross(base_vector_0, backbone_vector_0)
                aux_vector_b_0 = np.cross(aux_vector_a_0, base_vector_0)

                # return 3 orthogonal vectors in nucleotide, for quaternion
                vect_list_0 = (base_vector_0/np.linalg.norm(base_vector_0), \
                             aux_vector_a_0/np.linalg.norm(aux_vector_a_0), \
                             aux_vector_b_0/np.linalg.norm(aux_vector_b_0))

                for n, nucl in enumerate(self.oligos_list[o][s]):
                    [vh_1, index_1, is_fwd_1] = self.oligos_list_to_nucleotide_info( o, s, n)
                    [axis_1, backbone_1]      = self.nucleotide_matrix[vh_1][index_1][is_fwd_1].position

                    if n < len(self.oligos_list[o][s]) - 1:
                        [vh_2, index_2, is_fwd_2] = self.oligos_list_to_nucleotide_info(o, s, n + 1)
                        [axis_2, backbone_2]      = self.nucleotide_matrix[vh_2][index_2][is_fwd_2].position

                        base_vector_1     = axis_1     - backbone_1
                        backbone_vector_1 = backbone_2 - backbone_1
                    elif n == len(oligos_list[o][s]) - 1:
                        [vh_2, index_2, is_fwd_2] = self.oligos_list_to_nucleotide_info(o, s, n - 1)
                        [axis_2, backbone_2]      = self.nucleotide_matrix[vh_2][index_2][is_fwd_2].position
                        base_vector_1     = axis_1 - backbone_1
                        backbone_vector_1 = - (backbone_2 - backbone_1)

                    aux_vector_a_1 = np.cross(base_vector_1, backbone_vector_1)
                    aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
                    vect_list_1 = (base_vector_1+np.array([0.00001,0,0])/np.linalg.norm(base_vector_1+np.array([0.00001,0,0])), \
                                 aux_vector_a_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_a_1+np.array([0.00001,0,0])), \
                                 aux_vector_b_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_b_1+np.array([0.00001,0,0])))

                    nucl = self.nucleotide_matrix[vh_1][index_1][is_fwd_1]
                    nucl.vectors_body_frame  = vect_list_1
                    nucl.points_global_frame = [backbone_1, axis_1, aux_vector_a_1 + backbone_1]
                    nucl_quaternion          = vectortools.systemQuaternion(vect_list_0, vect_list_1)
                    nucl.quaternion          = [nucl_quaternion.w, \
                                                                              nucl_quaternion.x, \
                                                                              nucl_quaternion.y, \
                                                                              nucl_quaternion.z]

                    self.nucleotide_matrix[vh_1][index_1][is_fwd_1] = nucl      #This line is unnecessary

class DSNucleotide:
    '''
    Two sense/antisense nucleotide that are part of dsDNA strand
    '''
    def __init__(self):
        self.fwd_nucleotide     = None      # Nucleotide in forward direction (reference frame)
        self.rev_nucleotide     = None      # Nucleotide in reverse direction
        self.type               = 1
        self.visited            = False
        self.rigid              = False
        self.skip               = None

        #Connections
        self.rigid_connections  = []   #Rigid connections
        self.soft_connections   = []   #Soft connections

class SSNucleotide:
    '''
    Nucleotide that is part of ssDNA strand
    '''
    def __init__(self):
        self.nucleotide        = None
        self.type              = 0
        self.visited           = False
        self.rigid             = False
        self.skip              = None

        #Connections
        self.rigid_connections = []  #Rigid connections
        self.soft_connections  = []  #Soft connections

class Nucleotide:
    '''
    Fixed attributes of a nucleotide
    Attributes: index, position, strand, vh
    '''
    def __init__(self):
        self.direction                    = None                            # 1 is fwd, 0 is reverse
        self.is_fwd                       = None                            # 0: reverse, 1:forward
        self.index                        = None                            # z position in cadnano's unit
        self.strand                       = None                            # strand #
        self.vh                           = None                            # virtual helix this nucleotide belongs to
        self.skip                         = False                           # Skip value for the nucleotide

        self.points_global_frame          = None                            # backbone, sidechain and aux points in global frame
        self.quaternion                   = None                            # quaternion orientation for this nucleotide
        self.vectors_body_frame           = None                            # orthogonal vectors in the body reference frame for quaternion calculation
        self.position                     = None                            # Nucleotide position

        #Body/simulation variables
        self.body                         = None                            # body this nucleotide belongs to
        self.body_num                     = 0                               # body number
        self.simulation_nucleotide_num    = 0                               # body nucleotide number

class Body:
    '''
    Fixed attributes of a body.
    A body is a combination of neighboring vhs, making up a
    collection of nucleotides that move together
    Attributes: vhs, com_position, com_quaternion, nucleotides
    '''
    def __init__(self):
        self.comass_position   = None
        self.comass_quaternion = None
        self.moment_inertia    = None
        self.nucleotide_types  = []
        self.nucleotides       = []
        self.vhs               = []
        self.type              = None   #0: soft, 1:rigid

    def add_nucleotide_type(self, nucleotide_type):
        self.nucleotide_types.append(nucleotide_type)

    def add_vh(self, vh):
        self.vhs.add(vh)

    def initialize(self):
        '''
        Given a list of oligos, each composed of a list of strands
        each composed of a list of nucleotides,
        first populate each body's nucleotide and Vh and
        then calculate the other attributes.
        '''

        positions = [nucleotide.position[1] for nucleotide in self.nucleotides]
        self.comass_position    = vectortools.calculateCoM(positions)
        self.moment_inertia     = vectortools.calculateMomentInertia(positions)
        self.comass_quaternion  = [1., 0., 0., 0.]

class RigidBodySimulation:
    '''
    Rigid body simulation class for Cadnano designs
    '''

    def __init__(self):
        self.origami                 = None
        self.num_steps               = None
        self.ssDNA_harmonic_bond     = {'r0':None, 'k0':None}
        self.ssDNA_harmonic_angle    = {'a0':None, 'k0':None}

        self.dsDNA_harmonic_bond     = {'r0':None, 'k0':None}
        self.dsDNA_harmonic_angle    = {'a0':None, 'k0':None}

        self.bodies_comass_positions = []
        self.bodies_moment_inertia   = []

        self.snapshot                = None

        self.body_types              = []
        self.bond_types              = []

        #Rigid/soft bodies from Origami structure
        self.num_rigid_bodies        = 0
        self.num_soft_bodies         = 0
        self.rigid_bodies            = None
        self.soft_bodies             = None

    def initialize_relax_md(self):
        '''
        Initialize relaxation protocol
        '''
        context.initialize("");
        relax_sim = context.SimulationContext();

    def initialize_particles(self):
        '''
        Initialize particle positions, moment of inertia and velocities
        '''

        #Retrieve origami rigid body information
        self.rigid_bodies = self.origami.rigid_bodies
        self.soft_bodies  = self.origami.soft_bodies

        self.num_rigid_bodies = len(self.rigid_bodies)
        self.num_soft_bodies  = len(self.soft_bodies)

        self.rigid_bodies_comass_positions  = [body.comass_position for body in self.rigid_bodies]
        self.center_of_mass                 = np.average(np.asarray(self.rigid_bodies_comass_positions)[:,:3], axis=0)

        self.rigid_bodies_comass_positions -= self.center_of_mass
        self.rigid_bodies_moment_inertia    = [body.moment_inertia for body in self.rigid_bodies]

        self.soft_bodies_comass_positions  = [body.comass_position for body in self.soft_bodies]
        self.soft_bodies_comass_positions -= self.center_of_mass
        self.soft_bodies_moment_inertia    = [body.moment_inertia for body in self.soft_bodies]


        self.body_types  = ["rigid_body"+"_"+str(i) for i in range(self.num_rigid_bodies)]
        self.body_types += ["nucleotides"]

        self.snapshot = data.make_snapshot(N = self.num_rigid_bodies+self.num_soft_bodies,
                                          box = data.boxdim(Lx=120, Ly=120, Lz=300),
                                          particle_types = self.body_types,
                                          bond_types = ['interbody']);

        self.snapshot.particles.position[:]       = np.vstack((self.rigid_bodies_comass_positions,self.soft_bodies_comass_positions))
        self.snapshot.particles.moment_inertia[:] = np.vstack((self.rigid_bodies_moment_inertia  ,self.soft_bodies_moment_inertia))

        #particle types
        for i in range(self.num_rigid_bodies):
            self.snapshot.particles.typeid[i] = i

        #particle types
        for i in range(self.num_rigid_bodies,self.num_rigid_bodies+self.num_soft_bodies):
            self.snapshot.particles.typeid[i] = self.num_rigid_bodies

        self.snapshot.particles.velocity[:] = np.random.normal(0.0, np.sqrt(0.8 / 1.0), [self.snapshot.particles.N, 3]);

    def create_rigid_bodies(self):
        # Read the snapshot and create neighbor list
        self.system = init.read_snapshot(self.snapshot);
        self.nl     = md.nlist.stencil();

        # Create rigid particles
        self.rigid = md.constrain.rigid();
        for b, body in enumerate(self.rigid_bodies):
            body_type            = self.body_types[b]

            nucleotide_positions = [nucleotide.position[1] for nucleotide in body.nucleotides]
            #move particles to body reference frame
            nucleotide_positions -= body.comass_position
            self.rigid.set_param(body_type, \
                        types=['nucleotides']*len(nucleotide_positions), \
                        positions = nucleotide_positions);

        self.rigid.create_bodies()

    def create_bonds(self):
        '''
        Create interbody bonds
        '''
        self.nucleotide_bonds = self.origami.inter_nucleotide_connections

        for connection in self.nucleotide_bonds:
            delta = self.num_rigid_bodies
            nucleotide_num_1,nucleotide_num_2 = connection
            self.system.bonds.add('interbody', delta + nucleotide_num_1, delta + nucleotide_num_2)

    def set_harmonic_bonds(self):
        '''
        Set harmonic bonds
        '''
        self.harmonic = md.bond.harmonic()
        self.harmonic.bond_coeff.set('interbody', k=10.0    , r0=0.5);

        # fix diameters for vizualization
        for i in range(0, self.num_rigid_bodies):
            self.system.particles[i].diameter = 2.0
        for i in range(self.num_rigid_bodies, len(self.system.particles)):
            self.system.particles[i].diameter = 0.5

    def set_lj_potentials(self):
        '''
        Set LJ potentials
        '''
        wca = md.pair.lj(r_cut=2.0**(1/6), nlist=self.nl)
        wca.set_params(mode='shift')
        wca.pair_coeff.set(self.body_types, self.body_types, epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))

        ########## INTEGRATION ############
        md.integrate.mode_standard(dt=0.001, aniso=True);
        rigid     = group.rigid_center();
        non_rigid = group.nonrigid()
        combined  = group.union('combined',rigid,non_rigid)
        md.integrate.langevin(group=combined, kT=0.2, seed=42);

    def dump_settings(self,output_fname):
        '''
        Dump settings
        '''
        dump.gsd(output_fname,
                       period=1e4,
                       group=group.all(),
                       static=[],
                       overwrite=True);

    def update_quaternions(self):
        '''
        update global particle positions and quaternions
        '''
        for b, body in enumerate(self.rigid_bodies):
            for nucleotides in body.nucleotides:
                vh     = nucleotide.vh
                index  = nucleotide.index
                is_fwd = nucl.is_fwd
                simulation_num      = nucleotide.simulation_nucleotide_num
                nucleotide_position = self.system.particles[simulation_num].position

                nucl_quaternion_new = system.particles[simulation_num].orientation
                nucl_quaternion_new = vectortools.quat2Quat(nucl_quaternion_new)
                nucl_quaternion_old = nucleotide.quaternion
                nucl_quaternion_old = vectortools.quat2Quat(nucl_quaternion_old)

                quat = nucl_quaternion_new * nucl_quaternion_old
                quat = [quat.w, quat.x, quat.y, quat.z]
                self.origami.nucleotide_matrix[vh][index][is_fwd].position[1] = nucl_position
                self.origami.nucleotide_matrix[vh][index][is_fwd].quaternion  = quat

    def run(self,num_steps=1e6):
        run(num_steps)



def main():
    #Initialize cadnano
    app = cadnano.app()
    doc = app.document = Document()
    INPUT_FILENAME  = '../data/Hinge_v5.1-L1-R1.json'
    OUTPUT_FILENAME = '../data/Hinge_v5.1-L1-R1.gsd'

    doc.readFile(INPUT_FILENAME);

    #Parse the structure for simulation
    new_origami      = Origami()
    new_origami.part = doc.activePart()
    new_origami.list_oligos()
    new_origami.initialize_nucleotide_matrix()
    new_origami.find_skips()
    new_origami.create_oligos_list()
    new_origami.get_connections()
    new_origami.assign_nucleotide_types()
    new_origami.incorporate_skips()
    new_origami.assign_nucleotide_connections()
    new_origami.cluster_into_bodies()
    new_origami.parse_soft_connections()


    #Prepare the simulation
    new_simulation         = RigidBodySimulation()
    new_simulation.origami = new_origami
    new_simulation.initialize_relax_md()
    new_simulation.initialize_particles()
    new_simulation.create_rigid_bodies()
    new_simulation.create_bonds()
    new_simulation.set_harmonic_bonds()
    new_simulation.set_lj_potentials()
    new_simulation.dump_settings(OUTPUT_FILENAME)
    new_simulation.run(1e6)


if __name__ == "__main__":
  main()
