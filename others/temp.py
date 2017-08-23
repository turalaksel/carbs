def generateVectorsandQuaternions(oligos_array):
    '''
    Given an array of oligos, fills in nucleotide attributes
    '''
    nucl_list_list = listOflistsOfNucleotides(oligos_array)

    for o, oligo in enumerate(nucl_list_list):
        [base_0, bckbone_0] = nucl_list_list[0][0].position
        [base_1, bckbone_1] = nucl_list_list[0][1].position

        base_vector_0 = base_0 - bckbone_0
        polymer_vector_0 = bckbone_1 - bckbone_0
        aux_vector_a_0 = np.cross(base_vector_0, polymer_vector_0)
        aux_vector_b_0 = np.cross(aux_vector_a_0, base_vector_0)

        vect_list_0 = (base_vector_0/np.linalg.norm(base_vector_0), \
                     aux_vector_a_0/np.linalg.norm(aux_vector_a_0), \
                     aux_vector_b_0/np.linalg.norm(aux_vector_b_0))

        for n, nucl in enumerate(nucl_list_list[o]):
            [base_1, bckbone_1] = nucl_list_list[o][n].position
            if n < len(nucl_list_list[o]) - 1:
                [base_2, bckbone_2] = nucl_list_list[o][n + 1].position
                base_vector_1 = base_1 - bckbone_1
                polymer_vector_1 = base_2 - base_1
            elif n == len(nucl_list_list[o]) - 1:
                [base_2, bckbone_2] = nucl_list_list[o][n - 1].position
                base_vector_1 = base_1 - bckbone_1
                polymer_vector_1 = - (base_2 - base_1)

            aux_vector_a_1 = np.cross(base_vector_1, polymer_vector_1)
            aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
            vect_list_1 = (base_vector_1/np.linalg.norm(base_vector_1), \
                         aux_vector_a_1/np.linalg.norm(aux_vector_a_1), \
                         aux_vector_b_1/np.linalg.norm(aux_vector_b_1))

            nucl_quaternion = vTools.systemQuaternion(vect_list_0, vect_list_1)
            nucl.add_quaternion(nucl_quaternion.w)
            nucl.add_quaternion(nucl_quaternion.x)
            nucl.add_quaternion(nucl_quaternion.y)
            nucl.add_quaternion(nucl_quaternion.z)
            nucl.add_vectors(vect_list_1)

    return(nucl_list_list)
