#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2017-12-09 13:25:51
# @Author  : Tural Aksel (turalaksel@gmail.com)
# @Link    : http://example.org
# @Version : $Id$

import numpy as np
import cadnano
import vectortools
import functools
import argparse
import os

from cadnano.document import Document
from hoomd import *
from hoomd import md

import gsd.hoomd

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
        self.vh_vh_crossovers       = None   # Given vh1, vh2, vh_vh_crossovers[vh_1][vh_2] is the number of xovers between them
        self.long_range_connections = {}     # Dict connecting pointer_1 (vh, index, is_fwd) to pointer_2
        self.short_range_connections= {}     # Dict connecting pointer_1 (vh, index, is_fwd) to pointer_2
        self.soft_connections       = {}     # Dict of pointers referring to nucleotides separated by skip

        #Distance constraints
        self.crossover_distance      = 2.0   # Distance in nm

    def parse_soft_connections(self):
        self.inter_rigid_body_connections = set()
        self.inter_nucleotide_connections = set()

        for pointer_1, pointer_2 in self.soft_connections.items():
            vh_1, index_1, is_fwd_1 = pointer_1
            vh_2, index_2, is_fwd_2 = pointer_2

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
        Incorporate skips in the model by creating 'short_range_connections'
        between nucleotides in the beginning and end of skip region
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

        #2. Make short range connection between beginning and end of skip
        for vh, is_fwd, idx_begin, idx_end in self.skip_boundaries:
            pointer_1 = (vh, idx_begin - 1, 1)
            pointer_2 = (vh, idx_end + 1, 1)

            if not is_fwd:
                pointer_1 = (vh, idx_end + 1, 0)
                pointer_2 = (vh, idx_begin - 1, 0)

            self.short_range_connections[pointer_1] = pointer_2

    def assign_nucleotide_types(self):
        '''
        Build the nucleotide network for an origami design
        '''
        self.nucleotide_type_list = []
        for vh in range(len(self.nucleotide_matrix)):
            for idx in range(len(self.nucleotide_matrix[vh])):

                #Determine the nucleotide type (ssNucleotide vs dsNucletide) and nucleotide connections
                current_nucleotide_rev = self.nucleotide_matrix[vh][idx][0]
                current_nucleotide_fwd = self.nucleotide_matrix[vh][idx][1]

                if current_nucleotide_fwd != None and current_nucleotide_rev != None:
                    ds_nucleotide = DSNucleotide()
                    ds_nucleotide.fwd_nucleotide   = current_nucleotide_fwd
                    ds_nucleotide.rev_nucleotide   = current_nucleotide_rev
                    ds_nucleotide.type             = 1
                    ds_nucleotide.skip             = current_nucleotide_fwd.skip
                    self.nucleotide_type_matrix[vh][idx] = ds_nucleotide
                    self.nucleotide_type_list.append(ds_nucleotide)

                elif current_nucleotide_fwd != None:
                    ss_nucleotide                  = SSNucleotide()
                    ss_nucleotide.nucleotide       = current_nucleotide_fwd
                    ss_nucleotide.type             = 0
                    ss_nucleotide.skip             = current_nucleotide_fwd.skip

                    self.nucleotide_type_matrix[vh][idx] = ss_nucleotide
                    self.nucleotide_type_list.append(ss_nucleotide)

                elif current_nucleotide_rev != None:

                    ss_nucleotide                  = SSNucleotide()
                    ss_nucleotide.nucleotide       = current_nucleotide_rev
                    ss_nucleotide.type             = 0
                    ss_nucleotide.skip             = current_nucleotide_rev.skip

                    self.nucleotide_type_matrix[vh][idx] = ss_nucleotide
                    self.nucleotide_type_list.append(ss_nucleotide)

    def assign_nucleotide_connections(self):
        '''
        Assign nucleotide connections by looping over nucleotides in nucleotide_type_matrix
        then appending to nucleotide_type_matrix[vh][idx].rigid_connections
        the nucleotide_type_matrix[vh][idx] for both next and previous idx, if existent.
        '''
        #1. Add the base-stacking interactions
        for vh in range(len(self.nucleotide_type_matrix)):
            for idx in range(len(self.nucleotide_type_matrix[vh])):

                if self.nucleotide_type_matrix[vh][idx] == None or self.nucleotide_type_matrix[vh][idx].skip:
                    continue

                self.nucleotide_type_matrix[vh][idx].rigid_connections = []
                self.nucleotide_type_matrix[vh][idx].soft_connections  = []

                #Get the type for nucleotide (ssDNA or dsDNA?)
                type_1 = self.nucleotide_type_matrix[vh][idx].type

                #Pointer 1
                pointer1_rev = (vh, idx, 0)
                pointer1_fwd = (vh, idx, 1)

                # Calculate connections between idx and next idx, if existent
                if idx+1 < len(self.nucleotide_type_matrix[vh]) \
                   and not self.nucleotide_type_matrix[vh][idx+1] == None \
                   and not self.nucleotide_type_matrix[vh][idx+1].skip:

                    type_2 = self.nucleotide_type_matrix[vh][idx+1].type

                    if type_1*type_2: # both types are DSNucleotide, make the connection RIGID
                        self.nucleotide_type_matrix[vh][idx].rigid_connections.append(self.nucleotide_type_matrix[vh][idx+1])
                    else: # at least one is SSNucleotide, make a soft connection either in the fwd or rev direction
                        self.nucleotide_type_matrix[vh][idx].soft_connections.append(self.nucleotide_type_matrix[vh][idx+1])

                        #Pointer 2
                        pointer2_rev = (vh, idx+1, 0)
                        pointer2_fwd = (vh, idx+1, 1)

                        if pointer1_fwd in self.short_range_connections.keys() and self.short_range_connections[pointer1_fwd] == pointer2_fwd:
                            self.soft_connections[pointer1_fwd] = pointer2_fwd

                        elif self.short_range_connections[pointer2_rev] == pointer1_rev: #ssDNA connection is in the reverse direction
                            self.soft_connections[pointer2_rev] = pointer1_rev

                # Calculate connections between idx and previous idx, if existent
                if idx-1 >= 0 \
                   and not self.nucleotide_type_matrix[vh][idx-1] == None \
                   and not self.nucleotide_type_matrix[vh][idx-1].skip:

                    type_2 = self.nucleotide_type_matrix[vh][idx-1].type

                    if type_1*type_2:
                        self.nucleotide_type_matrix[vh][idx].rigid_connections.append(self.nucleotide_type_matrix[vh][idx-1])
                    else:
                        self.nucleotide_type_matrix[vh][idx].soft_connections.append(self.nucleotide_type_matrix[vh][idx-1])

                        pointer2_fwd = (vh, idx-1, 1)
                        pointer2_rev = (vh, idx-1, 0)

                        if pointer2_fwd in self.short_range_connections.keys() and self.short_range_connections[pointer2_fwd] == pointer1_fwd:
                            self.soft_connections[pointer2_fwd] = pointer1_fwd

                        elif self.short_range_connections[pointer1_rev] == pointer2_rev:
                            self.soft_connections[pointer1_rev] = pointer2_rev

        #2. Add short range connections that are not adjacent in sequence due to skips
        for pointer1, pointer2 in self.short_range_connections.items():
            vh1, idx1, is_fwd1 = pointer1
            vh2, idx2, is_fwd2 = pointer2

            #If the bases are not adjacent in sequence, add the connections to soft connections
            if abs(idx1-idx2) > 1:
                #Add the connections first in nucleotide type matrix
                self.nucleotide_type_matrix[vh1][idx1].soft_connections.append(self.nucleotide_type_matrix[vh2][idx2])
                self.nucleotide_type_matrix[vh1][idx1].soft_connections.append(self.nucleotide_type_matrix[vh2][idx2])
                self.soft_connections[pointer1] = pointer2

        #3. Add the crossover connections
        for pointer_1, pointer_2 in self.crossovers.items():
            (vh_1, index_1, is_fwd_1) = pointer_1
            (vh_2, index_2, is_fwd_2) = pointer_2

            type_1 = self.nucleotide_type_matrix[vh_1][index_1].type
            type_2 = self.nucleotide_type_matrix[vh_2][index_2].type

            if self.vh_vh_crossovers[vh_1][vh_2] > 1 and type_1*type_2: #make rigid if more than 1 xover and both nucleotides are dsDNA
                self.nucleotide_type_matrix[vh_1][index_1].rigid_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
                self.nucleotide_type_matrix[vh_2][index_2].rigid_connections.append(self.nucleotide_type_matrix[vh_1][index_1])
            else: # make soft otherwise
                self.nucleotide_type_matrix[vh_1][index_1].soft_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
                self.nucleotide_type_matrix[vh_2][index_2].soft_connections.append(self.nucleotide_type_matrix[vh_1][index_1])

                #Add the connection to soft connection list
                self.soft_connections[pointer_1] = pointer_2

        #4. Add long-range connections (always soft!)
        for pointer_1, pointer_2 in self.long_range_connections.items():
            (vh_1, index_1, is_fwd_1) = pointer_1
            (vh_2, index_2, is_fwd_2) = pointer_2

            self.nucleotide_type_matrix[vh_1][index_1].soft_connections.append(self.nucleotide_type_matrix[vh_2][index_2])
            self.nucleotide_type_matrix[vh_2][index_2].soft_connections.append(self.nucleotide_type_matrix[vh_1][index_1])

            #Add the connection to soft connection list
            self.soft_connections[pointer_1] = pointer_2

    def get_connections(self):
        '''
        Populate 3' connections for each (staple / scaffold) strand
        '''
        self.crossovers             = {}
        self.long_range_connections = {}
        for vh in range(self.num_vhs):
            staple_strandSet   = self.part.getStrandSets(vh)[not(vh % 2)]
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

        #1. Identify the clusters using depth-first-search
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
        Creates an empty matrix of len = vh_length x index_length
        to be populated with all nucleotides in part (fwd and rev)
        It has the form nucleotide_matrix[vh][index][rev or fwd]
        '''
        self.num_vhs = len(list(self.part.getIdNums()))
        num_bases    = self.part.getVirtualHelix(0).getSize()

        self.nucleotide_matrix       = [[[None,None]  for idx in range(num_bases)] for vh in range(self.num_vhs)]
        self.nucleotide_type_matrix  = [[ None  for idx in range(num_bases)] for vh in range(self.num_vhs)]
        self.vh_vh_crossovers        = [[0  for vh in range(self.num_vhs)] for vh in range(self.num_vhs)]
        self.skip_matrix             = [[[False,False]  for idx in range(num_bases)] for vh in range(self.num_vhs)]

    def connection3p(self, strand):
        '''
        Given a strand, returns the vhelix to which the 3p end
        connects to, if the distance is not too far
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

            if distance < self.crossover_distance:
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

    def create_strand_list_and_populate_nucleotide_matrix(self, oligo):
        '''
        Given an oligo, returns a list of strands,
        each containing the pointers ([vh][index][is_fwd]) to the
        nucleotides making up such strand and *populate nucleotides matrix*
        with attributes for this oligo
        '''
        if self.nucleotide_matrix == None:
            self.initialize_nucleotide_matrix()

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
        Returns a tuple list (aka pointer) [vh, index, is_fwd]
        for the nucleotide oligos_list[i][j][k] where
        i,j,k are the indices for the oligo, strand, and nucleotide.
        '''

        [vh, index, is_fwd] = self.oligos_list[i][j][k]
        return [vh, index, is_fwd]

    def determine_coordinate_limits(self):
        '''
        Determine origami coordinate limits 
        '''
        coordinates = []
        for nucleotide in self.nucleotide_list:
            coordinates.append(nucleotide.position[1])

        #Make the coordinates array
        coordinates = np.array(coordinates)

        #Determine min and max along all dimensions
        self.min_x  = np.min(coordinates[:,0])
        self.max_x  = np.max(coordinates[:,0])

        self.min_y  = np.min(coordinates[:,1])
        self.max_y  = np.max(coordinates[:,1])

        self.min_z  = np.min(coordinates[:,2])
        self.max_z  = np.max(coordinates[:,2])

    def create_oligos_list(self):
        '''
        Given an array of oligos in part, returns a list of oligos,
        each containing a list of strands, each containing a
        list of nucleotides making up the part.
        In the process, also populate the nucleotide_matrix w/ nucleotides
        '''
        self.oligos_list     = []
        self.nucleotide_list = []
        for oligo in self.oligos:
            strand_list = self.create_strand_list_and_populate_nucleotide_matrix(oligo)
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

class DSNucleotide:
    '''
    Fwd and Rev (sense/antisense) nucleotides making up a double strand nucleotide
    '''
    def __init__(self):
        self.fwd_nucleotide     = None                 # Nucleotide in forward direction (reference frame)
        self.rev_nucleotide     = None                 # Nucleotide in reverse direction
        self.type               = 1                    # 1: double strand (dsDNA)
        self.visited            = False                # To be used by depth-first-search
        self.rigid              = False
        self.skip               = None                 # Whether this nucleotide is a skip in cadnano

        #Connections
        self.rigid_connections  = []                   # List of rigid connections
        self.soft_connections   = []                   # List of soft connections

class SSNucleotide:
    '''
    Nucleotide that is part of ssDNA strand
    '''
    def __init__(self):
        self.nucleotide        = None
        self.type              = 0                    # 0: single strand (ssDNA)
        self.visited           = False                # To be used by depth-first-search
        self.rigid             = False
        self.skip              = None                 # Whether this nucleotide is a skip in cadnano

        #Connections
        self.rigid_connections = []                   # List of rigid connections
        self.soft_connections  = []                   # List of soft connections

class Nucleotide:
    '''
    Fixed attributes of a nucleotide
    '''
    def __init__(self):
        self.direction                    = None      # 1 is fwd, 0 is reverse
        self.is_fwd                       = None      # 0: reverse, 1:forward
        self.index                        = None      # z position in cadnano's unit
        self.strand                       = None      # Nucleotide's strand number
        self.vh                           = None      # Nucleotide's virtual helix
        self.skip                         = False     # Skip value for the nucleotide
        self.position                     = None      # Nucleotide position

        # Body / simulation variables
        self.body                         = None      # body (class) this nucleotide belongs to
        self.body_num                     = 0         # body number
        self.simulation_nucleotide_num    = 0         # nucleotide number wrt the hoomd simulation

class Body:
    '''
    Fixed attributes of a body.
    A body is a combination of neighboring vhs that move together during
    relaxation as one rigid body. HOOMD will need its center of mass position and moment of inertia.
    '''
    def __init__(self):
        self.comass_position   = None                 # position of body's center of mass
        self.moment_inertia    = None                 # body's moment of intertia (calculated via vectortools)
        self.nucleotide_types  = []                   # list of nucleotide types belonging to this body
        self.nucleotides       = []                   # list of nucleotides belonging to this body
        self.vhs               = []                   # list of vhs belonging to this body
        self.type              = None                 # 0: soft, 1:rigid

    def add_nucleotide_type(self, nucleotide_type):
        self.nucleotide_types.append(nucleotide_type)

    def add_vh(self, vh):
        self.vhs.add(vh)

    def initialize(self):
        '''
        Given a collection of nucleotides making up a body, initialize the body
        by calculating the following properties: comass_position, and moment_inertia
        '''

        # extract the position of backbone bead (1) acquired from cadnano
        positions = [nucleotide.position[1] for nucleotide in self.nucleotides]
        self.comass_position    = vectortools.calculate_comass(positions)
        self.moment_inertia     = vectortools.calculate_moment_inertia(positions)

class RigidBodySimulation:
    '''
    Rigid body simulation class for Cadnano designs
    '''

    def __init__(self):
        self.origami                  = None
        self.num_steps                = None
        self.ssDNA_harmonic_bond      = {'r0':None, 'k0':None}
        self.ssDNA_harmonic_angle     = {'a0':None, 'k0':None}

        self.dsDNA_harmonic_bond      = {'r0':None, 'k0':None}
        self.dsDNA_harmonic_angle     = {'a0':None, 'k0':None}

        self.interbody_harmonic_bond  = {'r0':None, 'k0':None}
        self.interbody_harmonic_angle = {'a0':None, 'k0':None}

        self.bodies_comass_positions  = []
        self.bodies_moment_inertia    = []

        self.snapshot                 = None

        self.body_types               = []
        self.bond_types               = []

        #Simulation constant for determining spring constant and time step
        #Multiplication of maximum_bond_distance*spring_constant*(time_step**2) should be smaller than this value
        self.SIMULATION_CONSTANT      = 0.000022

        self.simulation_step_cutoff   = 1E-4

        #Rigid/soft bodies from Origami structure
        self.num_rigid_bodies         = 0
        self.num_soft_bodies          = 0
        self.rigid_bodies             = None
        self.soft_bodies              = None

    def determine_box_dimensions(self,scale=6):
        '''
        Estimate the box dimensions from particle coordinates
        '''
        self.Lx = scale*int((self.origami.max_x - self.origami.min_x)/2.0)
        self.Ly = scale*int((self.origami.max_y - self.origami.min_y)/2.0)
        self.Lz = scale*int((self.origami.max_z - self.origami.min_z)/2.0)

        print('Box dimensions:%d-%d-%d'%(self.Lx,self.Ly,self.Lz))

    def set_interbody_harmonic_bond(self,r0=0.745,k0=1.0):
        '''
        Set interbody harmonic bond parameters
        '''
        self.interbody_harmonic_bond['r0'] = r0
        self.interbody_harmonic_bond['k0'] = k0

    def initialize_relax_md(self):
        '''
        Initialize relaxation protocol
        '''
        context.initialize("");
        relax_sim = context.SimulationContext();

    def read_gsd(self):
        '''
        Read gsd file
        '''
        self.trajectory     = gsd.hoomd.open(self.gsd_filename,'rb')
        self.final_snapshot = self.trajectory[-1]

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

        self.snapshot = data.make_snapshot(N = self.num_rigid_bodies + self.num_soft_bodies,
                                          box = data.boxdim(Lx=self.Lx, Ly=self.Ly, Lz=self.Lz),
                                          particle_types = self.body_types,
                                          bond_types = ['interbody']);

        self.snapshot.particles.position[:]       = np.vstack((self.rigid_bodies_comass_positions, self.soft_bodies_comass_positions))
        self.snapshot.particles.moment_inertia[:] = np.vstack((self.rigid_bodies_moment_inertia  , self.soft_bodies_moment_inertia))

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
            
            #Move particles to body reference frame
            nucleotide_positions -= body.comass_position
            self.rigid.set_param(body_type,
                                types=['nucleotides']*len(nucleotide_positions),
                                positions = nucleotide_positions);

        self.rigid.create_bodies()

    def take_snapshot(self):
        '''
        Take system synaphot
        '''
        self.final_snapshot = self.system.take_snapshot()

    def create_bonds(self):
        '''
        Create interbody bonds
        '''
        self.nucleotide_bonds = self.origami.inter_nucleotide_connections

        #Nucleotide number offset
        delta = self.num_rigid_bodies
        for connection in self.nucleotide_bonds:
            nucleotide_num_1, nucleotide_num_2 = connection
            self.system.bonds.add('interbody', delta + nucleotide_num_1, delta + nucleotide_num_2)

    def set_harmonic_bonds(self,spring_constant=1.0, spring_distance=0.745):
        '''
        Set harmonic bonds
        '''
        self.harmonic = md.bond.harmonic()
        self.harmonic.bond_coeff.set('interbody', k=spring_constant , r0=spring_distance);

        #Save the interbody harmonic bond parameters
        self.set_interbody_harmonic_bond(spring_distance,spring_constant)
        
        #Fix diameters for visualization
        for i in range(0, self.num_rigid_bodies):
            self.system.particles[i].diameter = 2.0
        for i in range(self.num_rigid_bodies, len(self.system.particles)):
            self.system.particles[i].diameter = 0.1

    def set_simulation_step(self,k0=1.0):
        '''
        Set simulation step for a spring constant
        '''
        self.simulation_step = np.sqrt(self.SIMULATION_CONSTANT/(k0*(self.max_bond_distance-self.interbody_harmonic_bond['r0'])))
        self.spring_constant = k0
        self.harmonic.bond_coeff.set('interbody', k=k0);
        md.integrate.mode_standard(dt=self.simulation_step, aniso=True);


    def get_distances(self):
        '''
        Get the interbody distanes
        '''
        #Nucleotide number offset
        delta = self.num_rigid_bodies

        #Initialize bond distances
        self.bond_distances = []

        for connection in self.nucleotide_bonds:
            #Get the nucleotide numbers
            nucleotide_num_1, nucleotide_num_2 = connection

            #Get the particle numbers in the simulation
            particle_num_1 = nucleotide_num_1 + delta
            particle_num_2 = nucleotide_num_2 + delta

            #Measure the distance
            particle_coor_1 = self.final_snapshot.particles.position[particle_num_1]
            particle_coor_2 = self.final_snapshot.particles.position[particle_num_2]
            
            distance = vectortools.eular_distance(particle_coor_1,particle_coor_2)

            self.bond_distances.append([particle_num_1,particle_num_2,distance])
        
        #Make the bond distances array
        self.bond_distances    = np.array(self.bond_distances)

        #Get the maximum value
        self.max_bond_distance = np.max(self.bond_distances[:,2])
        self.min_bond_distance = np.min(self.bond_distances[:,2])
        
    def set_lj_potentials(self):
        '''
        Set LJ potentials
        '''
        wca = md.pair.lj(r_cut=2.0**(1/6), nlist=self.nl)
        wca.set_params(mode='shift')
        wca.pair_coeff.set(self.body_types, self.body_types, epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))

    def set_simulation_settings(self,time_step=0.001, kT=1.0, rand_seed=42): 
        md.integrate.mode_standard(dt=time_step, aniso=True);
        rigid     = group.rigid_center();
        non_rigid = group.nonrigid()
        combined  = group.union('combined',rigid, non_rigid)
        md.integrate.langevin(group=combined, kT=kT, seed=rand_seed);

    def dump_settings(self,output_fname,save_period=1e3):
        '''
        Dump settings
        '''
        #Assign the output filename
        self.gsd_filename = output_fname 
        
        dump.gsd(output_fname,
                       period=save_period,
                       group=group.all(),
                       static=[],
                       overwrite=True);

    def run(self,num_steps=1e6):
        run(num_steps)

    def run_up_to(self,num_steps=1e4):
        run_upto(num_steps)

    def folding_protocol(self,simulation_time=1e4,num_itr=10):
        '''
        Folding protocol
        '''

        print('Starting maximum distance:%f'%(self.max_bond_distance))
        
        for itr in range(num_itr):
            
            #Ramping up the spring constant
            start_1  = itr*10
            finish_1 = (itr+1)*10
            for i in range(start_1,finish_1):
                self.run_up_to(simulation_time*(i+1))
                self.read_gsd()
                self.get_distances()
                self.set_simulation_step(k0=10**(itr)+(i-start_1)*10**itr)
                
                print('Folding step:%d, Maximum Distance:%f, Time Step:%e, Spring Constant:%f'%(i,self.max_bond_distance, self.simulation_step, self.spring_constant))
                
                #If simulation step is below the cutoff value, exit
                if self.simulation_step < self.simulation_step_cutoff:
                    sys.exit('Simulation step is below cutoff value!')    

            #Spring constant is kept constant
            start_2  = (itr+1)*10
            finish_2 = (itr+2)*10
            for i in range(start_2,finish_2):
                self.run_up_to(simulation_time*(i+1))
                self.read_gsd()
                self.get_distances()
                self.set_simulation_step(k0=10**(itr+1))
                
                print('Folding step:%d, Maximum Distance:%f, Time Step:%e, Spring Constant:%f'%(i,self.max_bond_distance, self.simulation_step, self.spring_constant))

                #If simulation step is below the cutoff value, exit
                if self.simulation_step < self.simulation_step_cutoff:
                    sys.exit('Simulation step is below cutoff value!')
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",  type=str,   help="Cadnano json file" )
    parser.add_argument("-o", "--output", type=str,   help="Output directory",          default='.')

    args = parser.parse_args()
    
    #Assign the parameters
    input_filename   = args.input
    output_directory = args.output

    #Check if input file exists
    if not os.path.isfile(input_filename):
        sys.exit('Input file does not exist!')

    #Check if output directory exists
    if not os.path.isdir(output_directory):
        sys.exit('Output directory does not exist!')

    #Define output filename
    head, tail      = os.path.split(input_filename)
    root, ext       = os.path.splitext(tail)
    output_filename = output_directory+'/'+root+'.gsd'

    #Initialize cadnano
    app = cadnano.app()
    doc = app.document = Document()

    #Read cadnano input file
    doc.readFile(input_filename);

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
    new_origami.determine_coordinate_limits()

    #Prepare the simulation
    new_simulation         = RigidBodySimulation()
    new_simulation.origami = new_origami
    new_simulation.determine_box_dimensions()
    new_simulation.initialize_relax_md()
    new_simulation.initialize_particles()
    new_simulation.create_rigid_bodies()
    new_simulation.create_bonds()
    new_simulation.set_harmonic_bonds()
    new_simulation.set_lj_potentials()
    new_simulation.set_simulation_settings()
    new_simulation.dump_settings(output_filename)
    new_simulation.take_snapshot()
    new_simulation.get_distances()
    new_simulation.set_simulation_step()
    new_simulation.folding_protocol()
    
if __name__ == "__main__":
  main()
