'''
This program returns a list of particle
positions and types from a PDB file

For that end, we have created Beads that represent
the atoms, Nucleotides, made out of beads, and Chains,
made out of nucleotides.
'''

import gsd
import gsd.hoomd
import gsd.pygsd
import numpy as np

class Bead:
    """
    Fixed attributes of a spherical particle
    Attributes: diameter, position, type.

    """
    def __init__(self):
        self.body      = []
        self.diameter  = []
        self.position  = []
        self.type	   = []
    def add_body(self, bead_body):
        self.body.append(bead_body)
    def add_diameter(self, bead_diameter):
        self.diameter.append(bead_diameter)
    def add_position(self, bead_position):
        self.position.append(bead_position)
    def add_type(self, bead_type):
        self.type.append(bead_type)

class Nucleotide:
    """
    Fixed attributes of a group of spherical particles making a nucleotide.
    Attributes: beads.

    """
    def __init__(self):
        self.beads = []
    def add_beads(self, nucleotide_beads):
        self.beads.append(nucleotide_beads)

class Chain:
    """
    Fixed attributes of a group of nucleotides tied to each other.
    Attributes: nucleotides.

    """
    def __init__(self):
        self.nucleotides = []
    def add_nucleotides(self, chain_nucleotides):
        self.nucleotides.append(chain_nucleotides)

def nucleotideCofM(nucleotide):
    center_of_mass = np.array([0.,0.,0.])
    positions_list = [b.position[0] for b in nucleotide.beads]
    positions_list_array = np.asarray(positions_list)
    # change here if particles have different masses !
    center_of_mass = np.average(positions_list_array[:,:3], axis=0)
    return(list(center_of_mass))

def shift2Center(list_of_chains):
    center = np.average(np.asarray(array)[:,:3], axis=0)
    array = np.asarray(array) - center
    return(array)

def getChainsFromPDB(input_file):
    previous_nucleotide_number = 0 # I assume the bodies are in order (1, 2...)
    list_of_chains = []
    new_chain = Chain()
    new_nucleotide = Nucleotide()

    for line in open(input_file):
        line = line.split()
        line_type = line[0]
        if line_type == 'ATOM': # we are still in the same chain
            body = int(line[5])
            if body != previous_nucleotide_number: #new nucleotide found
                new_chain.add_nucleotides(new_nucleotide)
                new_nucleotide = Nucleotide()
                previous_nucleotide_number += 1
            new_bead = Bead()
            new_bead.add_body(int(line[5]) - 1)
            new_bead.add_diameter(1.0)
            new_bead.add_position([float(line[6]),float(line[7]),float(line[8])])
            new_nucleotide.add_beads(new_bead)

        elif line_type == 'TER':
            list_of_chains.append(new_chain)
            new_chain = Chain()
    return(list_of_chains)
