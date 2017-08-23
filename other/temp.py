chains = fromPDB.getChainsFromPDB('doubleStrand.pdb')
for chain in chains:
    for nucleotide in chain.nucleotides:
        for bead in nucleotide.beads:
            print(bead.position[0])

bonds = [[n, min(n+1, 10)] for n in range(0, number_of_nucleotides - 1, 1)]


import numpy
