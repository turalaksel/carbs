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

def calculateCoM(positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the center of mass (vector3) of the system of particles
    Assumes equal masses
    '''
    center = np.average(np.asarray(positions)[:,:3], axis=0)
    return(center)

def calculateMomentInertia(positions):
    '''
    Given a list of arrays containing particle positions (vector3)
    Return the moment of inertia (vector3) of the system of particles
    Assumes equal masses
    '''
    inertia = np.array([0., 0., 0.])
    center = calculateCoM(positions)
    for p, pos in enumerate(positions):
        shifted_pos = pos - center
        new_inertia = np.multiply(shifted_pos, shifted_pos)
        inertia = inertia + new_inertia
    #re-scale particle masses so that body is not hugely slow
    #this needs to be tested
    inertia /= p
    return(inertia)

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
        body.com_position = calculateCoM(positions)
        body.moment_inertia = calculateMomentInertia(positions)
        body.com_quaternion = [1., 0., 0., 0.]
    return(bodies)

all_bodies = populateBody(list_of_list_of_nucleotides)

all_bodies[0].moment_inertia

vhs_of_body_list
nucl = list_of_list_of_nucleotides[1][0][0]
nucl.vh
for i in range(len(vhs_of_body_list)):
    if nucl.vh in vhs_of_body[i]:
        body_id = i



for vh in vhs_of_body_list: #vhs_of_body is a set!
    scaffold_strandSet = part.getStrandSets(vh)[vh % 2]
    staple_strandSet = part.getStrandSets(vh)[not(vh % 2)]

all_bodies = separateOrigamiParts(part)
list(all_bodies[0])


list_of_list_of_nucleotides = generateVectorsandQuaternions(oligos_array)

list_of_list_of_nucleotides[0][0][0].vh

bodies_list = [Body() for i in range(3)]
bodies_list[0].vhs = 0
bodies_list[0].vhs

if 0 in vhs_of_body[0]:
    print("yay")
