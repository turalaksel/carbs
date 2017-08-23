index_1st_nucl_in_strand = 0
for chain in range(len(list_of_list_of_nucleotides)):
    for strand in range(len(list_of_list_of_nucleotides[chain])):
        for nucl in range(len(list_of_list_of_nucleotides[chain][strand]) - 2):

            bckbone_1 = index_1st_nucl_in_strand + nucl
            sidechain_1 = total_num_nucl + 2*index_1st_nucl_in_strand + 2*nucl
            aux_1 = sidechain_1 + 1

            bckbone_2 = bckbone_1 + 1
            sidechain_2 = sidechain_1 + 1
            aux_2 = aux_1 + 1

            system.dihedrals.add('dihedral1',  sidechain_1, bckbone_1, bckbone_2, sidechain_2)
            system.dihedrals.add('dihedral21', sidechain_1, bckbone_1, aux_1, bckbone_2)
            system.dihedrals.add('dihedral22', sidechain_2, bckbone_2, aux_2, bckbone_1)
            system.dihedrals.add('dihedral31', aux_1, bckbone_1, sidechain_1, bckbone_2)
            system.dihedrals.add('dihedral32', aux_2, bckbone_2, sidechain_2, bckbone_1)


            print(chain, strand, nucl, bckbone_1, sidechain_1, aux_1)
        index_1st_nucl_in_strand += len(list_of_list_of_nucleotides[chain][strand])
