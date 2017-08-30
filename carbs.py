from os import path
import numpy as np
import cadnano
from cadnano.document import Document
import vector_tools as vTools #needed for vector / quaternion calculations

####### USER DEFINED AND GLOBAL VARIABLES ######
INPUT_FILENAME = 'input/tripod.json'
RELAX = False

global_nucl_matrix = []
global_connections_pairs = []

########### helper #############
# Helper functions for cadnano #
################################

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

def getNucleotide(nucl_pointers):
    '''
    Given a tuple of pointers in the form [vh, index, is_fwd],
    Returns the global nucleotide referent to the pointers
    '''
    [vh, index, is_fwd] = nucl_pointers
    return(global_nucl_matrix[vh][index][is_fwd])

def createOligosList(oligos_array, active_part):
    '''
    Given an array of oligos in part, returns a list of *oligos*,
    each containing a list of *strands(), each containing a
    list of *nucleotides() making up the part.
    '''
    global global_nucl_matrix
    oligos_list = []
    for oligo in oligos_array:
        strand_list = populateBasicNucleotideAttributes(oligo, active_part)
        oligos_list.append(strand_list)
    return(oligos_list)

def translateIndices(oligos_list, i, j, k):
    '''
    Oligos helper list gives a tuple list of pointers
    of the type [vh, index, is_fwd] which are needed to
    reference nucleotides in global_nucl_matrix
    this function translates oligo_helper_list into
    indices to be used by global_nucl_matrix
    '''
    [vh, index, is_fwd] = oligos_list[i][j][k]
    return([vh, index, is_fwd])

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

def populateAllNucleotideAttributes(oligos_list, active_part):
    '''
    Given an array of oligos, fills the remaining nucleotide attributes:
    vectors, quaternion, global_pts
    '''
    global global_nucl_matrix

    for o, oligo in enumerate(oligos_list):
        for s, strand in enumerate(oligo):

            [vh_0, index_0, is_fwd_0] = translateIndices(oligos_list, 0, 0, 0)
            [vh_1, index_1, is_fwd_1] = translateIndices(oligos_list, 0, 0, 1)

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

            for n, nucl in enumerate(oligos_list[o][s]):
                [vh_1, index_1, is_fwd_1] = translateIndices(oligos_list, o, s, n)
                [axis_1, backbone_1] = global_nucl_matrix[vh_1][index_1][is_fwd_1].position

                if n < len(oligos_list[o][s]) - 1:
                    [vh_2, index_2, is_fwd_2] = translateIndices(oligos_list, o, s, n + 1)
                    [axis_2, backbone_2] = global_nucl_matrix[vh_2][index_2][is_fwd_2].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = backbone_2 - backbone_1
                elif n == len(oligos_list[o][s]) - 1:
                    [vh_2, index_2, is_fwd_2] = translateIndices(oligos_list, o, s, n - 1)
                    [axis_2, backbone_2] = global_nucl_matrix[vh_2][index_2][is_fwd_2].position
                    base_vector_1 = axis_1 - backbone_1
                    backbone_vector_1 = - (backbone_2 - backbone_1)

                aux_vector_a_1 = np.cross(base_vector_1, backbone_vector_1)
                aux_vector_b_1 = np.cross(aux_vector_a_1, base_vector_1)
                vect_list_1 = (base_vector_1+np.array([0.00001,0,0])/np.linalg.norm(base_vector_1+np.array([0.00001,0,0])), \
                             aux_vector_a_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_a_1+np.array([0.00001,0,0])), \
                             aux_vector_b_1+np.array([0.00001,0,0])/np.linalg.norm(aux_vector_b_1+np.array([0.00001,0,0])))

                nucl = global_nucl_matrix[vh_1][index_1][is_fwd_1]
                nucl.vectors = vect_list_1
                nucl.global_pts = [backbone_1, axis_1, aux_vector_a_1 + backbone_1]
                nucl_quaternion = vTools.systemQuaternion(vect_list_0, vect_list_1)
                nucl.quaternion = [nucl_quaternion.w, \
                                                                          nucl_quaternion.x, \
                                                                          nucl_quaternion.y, \
                                                                          nucl_quaternion.z]
                global_nucl_matrix[vh_1][index_1][is_fwd_1] = nucl

############### relax ################
## Functions needed for relaxation  ##
######################################

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
            vh_1 = strand.idNum()
            index_1 = strand.idx3Prime()
            vh_2 = strand.connection3p().idNum()
            index_2 = strand.connection3p().idx5Prime()
            distance = distanceBetweenVhs(vh_1, index_1, vh_2, index_2)
            if distance < 6.0:
                return(vh_2)
            else:
                is_fwd_1 = int(strand.isForward())
                is_fwd_2 = int(strand.connection3p().isForward())
                conn_pointer_1 = [vh_1, index_1, is_fwd_1]
                conn_pointer_2 = [vh_2, index_2, is_fwd_2]
                if [conn_pointer_2, conn_pointer_1] not in global_connections_pairs:
                    global_connections_pairs.append([conn_pointer_1, conn_pointer_2])

def connection5p(strand):
    '''
    Given a strand, returns the vhelix to which the 5p end
    connects to, if the distance is not too far
    '''
    global global_connections_pairs
    if strand.connection5p() != None:
            vh_1 = strand.idNum()
            index_1 = strand.idx5Prime()
            vh_2 = strand.connection5p().idNum()
            index_2 = strand.connection5p().idx3Prime()
            distance = distanceBetweenVhs(vh_1, index_1, vh_2, index_2)
            if distance < 6.0:
                return(vh_2)
            else:
                is_fwd_1 = int(strand.isForward())
                is_fwd_2 = int(strand.connection5p().isForward())
                conn_pointer_1 = [vh_1, index_1, is_fwd_1]
                conn_pointer_2 = [vh_2, index_2, is_fwd_2]
                if [conn_pointer_2, conn_pointer_1] not in global_connections_pairs:
                    global_connections_pairs.append([conn_pointer_1, conn_pointer_2])

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

def populateBodiesNuclAndVhs(oligos_list):
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

    for o, oligo in enumerate(oligos_list):
        for s, strand in enumerate(oligo):
            for n, nucl in enumerate(strand):
                # seach in each body for this nucleotides' vh, assign to body
                [vh_1, index_1, is_fwd_1] = translateIndices(oligos_list, o, s, n)
                this_nucl = global_nucl_matrix[vh_1][index_1][is_fwd_1]
                body_id = [i for i in range(num_bodies) if this_nucl.vh in vhs_of_body[i]][0]
                this_nucl.body = body_id
                bodies[body_id].add_nucleotide(this_nucl)
                bodies[body_id].add_vh(this_nucl.vh)
    return(bodies)

def populateBody(oligos_list):
    '''
    Given a list of oligos, each composed of a list of strands
    each composed of a list of nucleotides,
    first populate each body's nucleotide and Vh and
    then calculate the other attributes.
    '''
    bodies = populateBodiesNuclAndVhs(oligos_list)

    for body in bodies:
        positions = [nucl.position[1] for nucl in body.nucleotides]
        body.com_position = vTools.calculateCoM(positions)
        body.moment_inertia = vTools.calculateMomentInertia(positions)
        body.com_quaternion = [1., 0., 0., 0.]
    return(bodies)

def findNuclPosition(nucl_pointer, bodies_list, global_ref = True):
    '''
    Given a nucleotide pointer = [vh, index, is_fwd]
    and a list of bodies already computed
    Returns the location of the nucleotide in the body
    if global = True, returns the location w.r.t.
    first nucleotide in body = 0
    '''
    my_nucl = getNucleotide(nucl_pointer)
    my_body = my_nucl.body
    position = [n for n, nucl in enumerate(bodies_list[my_body].nucleotides) \
                if nucl == my_nucl][0]
    if global_ref:
        for b in range(my_body):
            position += len(bodies_list[b].nucleotides)

    return(position)

############# cad #################
# start cadnano and read input file
###################################

app = cadnano.app()
doc = app.document = Document()
doc.readFile(INPUT_FILENAME);
part = doc.activePart()
oligos_array = oligosArray(part)

oligos_list = createOligosList(oligos_array, part)
populateAllNucleotideAttributes(oligos_list, part)
bodies = populateBody(oligos_list)

############# hoomd ################################################
# Start HOOMD code
####################################################################

from hoomd import *
from hoomd import md

# Start HOOMD in two different contexts
context.initialize("");
relax_sim = context.SimulationContext();
MD_sim = context.SimulationContext();

if RELAX == True:
    with relax_sim:
        num_rigid_bodies = len(bodies)

        bodies_com_positions = [body.com_position for body in bodies]
        bodies_com_positions -= np.average(np.asarray(bodies_com_positions)[:,:3], axis=0)

        bodies_mom_inertia = [body.moment_inertia for body in bodies]
        body_types = ["body"+"_"+str(i) for i in range(num_rigid_bodies)]
        body_types.append('nucleotides')

        snapshot = data.make_snapshot(N = num_rigid_bodies,
                                      box = data.boxdim(Lx=120, Ly=120, Lz=300),
                                      particle_types=body_types,
                                      bond_types = ['body', 'interbody']);

        snapshot.particles.position[:] = bodies_com_positions
        snapshot.particles.moment_inertia[:] = bodies_mom_inertia
        #particle types
        for i in range(len(bodies)):
            snapshot.particles.typeid[i] = i

        snapshot.particles.velocity[:] = np.random.normal(0.0,
                                         np.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

        # Bonds between connected rigid bodies
        def bodyConnections():
            connections = []
            for conn in global_connections_pairs:
                nucl0 = getNucleotide(conn[0])
                nucl1 = getNucleotide(conn[1])
                if [nucl1.body, nucl0.body] and [nucl0.body, nucl1.body] not in connections:
                    connections.append([nucl0.body, nucl1.body])
            return(connections)

        #create very weak bonds between com that are connected
        bonds = bodyConnections()
        snapshot.bonds.resize(len(bonds))
        snapshot.bonds.group[:] = bonds

        # Read the snapshot and create neighbor list
        system = init.read_snapshot(snapshot);
        nl = md.nlist.stencil();

        # Create rigid particles
        rigid = md.constrain.rigid();
        for b, body in enumerate(bodies):
            body_type = body_types[b]
            nucl_positions = [nucl.position[1] for nucl in body.nucleotides]
            # move particles to body reference frame
            nucl_positions -= body.com_position
            rigid.set_param(body_type, \
                        types=['nucleotides']*len(nucl_positions), \
                        positions = nucl_positions); #magic numbers. Check !!!

        rigid.create_bodies()

        harmonic = md.bond.harmonic()
        harmonic.bond_coeff.set('body', k=0.001, r0=10);
        harmonic.bond_coeff.set('interbody', k=0.1, r0=1.0)

        # fix diameters for vizualization
        for i in range(0, num_rigid_bodies):
            system.particles[i].diameter = 4.0
        for i in range(num_rigid_bodies, len(system.particles)):
            system.particles[i].diameter = 0.3

        ## Interbody bonds
        for conn in global_connections_pairs:
            nucl_0_body_location = findNuclPosition(conn[0], bodies, True)
            nucl_1_body_location = findNuclPosition(conn[1], bodies, True)
            delta = num_rigid_bodies

            system.bonds.add('interbody', delta + nucl_0_body_location, delta + nucl_1_body_location)

        ########## INTERACTIONS ############
        # LJ interactions
        wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
        wca.set_params(mode='shift')
        wca.pair_coeff.set(body_types, body_types, epsilon=1.0, sigma=1.0, r_cut=1.0*2**(1/6))

        ########## INTEGRATION ############
        md.integrate.mode_standard(dt=0.005, aniso=True);
        rigid = group.rigid_center();
        md.integrate.langevin(group=rigid, kT=0.5, seed=42);

        ########## DUMP & RUN ############
        dump.gsd("output/tripod_relax.gsd",
                       period=1e5,
                       group=group.all(),
                       static=[],
                       overwrite=True);

        run(2e5)

        #update global particle positions and quaternions
        i = num_rigid_bodies
        for b, body in enumerate(bodies):
            for nucl in body.nucleotides:
                vh = nucl.vh
                index = nucl.index
                is_fwd = int((nucl.direction + 1)/2)
                nucl_position = system.particles[i].position

                nucl_quaternion_new = system.particles[i].orientation
                nucl_quaternion_new = vTools.quat2Quat(nucl_quaternion_new)
                nucl_quaternion_old = nucl.quaternion
                nucl_quaternion_old = vTools.quat2Quat(nucl_quaternion_old)
                quat = nucl_quaternion_new * nucl_quaternion_old
                quat = [quat.w, quat.x, quat.y, quat.z]
                global_nucl_matrix[vh][index][is_fwd].position[1] = nucl_position
                global_nucl_matrix[vh][index][is_fwd].quaternion = quat
                i += 1


        import pickle
        with open('global.pckl', 'wb') as f:
            pickle.dump(global_nucl_matrix, f)

if RELAX == False:
    with MD_sim:
        #read global_nucl_matrix from pickle file
        import pickle
        f = open('global.pckl', 'rb')
        global_nucl_matrix = pickle.load(f)
        f.close()

        num_oligos = len(oligos_list)
        nucl_positions = [getNucleotide(pointer).position[1] for chain in oligos_list \
                            for strand in chain for pointer in strand]
        total_num_nucl = len(nucl_positions)

        snapshot = data.make_snapshot(N = total_num_nucl,
                                      box = data.boxdim(Lx=120, Ly=120, Lz=120),
                                      particle_types=['backbone','sidechain','aux'],
                                      bond_types = ['backbone','aux_sidechain'],
                                      dihedral_types = ['dihedral1', \
                                                        'dihedral21',\
                                                        'dihedral22',\
                                                        'dihedral31',\
                                                        'dihedral32']);

        # particle positions, types and moments of inertia
        nucl_positions = [getNucleotide(pointer).position[1] for chain in oligos_list \
                            for strand in chain for pointer in strand]
        nucl_quaternions = [getNucleotide(pointer).quaternion for chain in oligos_list \
                            for strand in chain for pointer in strand]

        snapshot.particles.position[:] = nucl_positions
        snapshot.particles.orientation[:] = nucl_quaternions
        snapshot.particles.moment_inertia[:] = [[1.,1.,1.]] #not correct. fix it
        snapshot.particles.typeid[:] = [0];

        # Backbone bonds
        bonds = []
        i = 0
        for chain in oligos_list:
            flat_chain = np.concatenate(chain)
            for n in range(i, i + len(flat_chain) - 1):
                bonds.append([n, n+1])
            i += len(flat_chain)

        #fix: add extra bond for circular staples / scaffold
        snapshot.bonds.resize(len(bonds))
        snapshot.bonds.group[:] = bonds

        # Read the snapshot and create neighbor list
        system = init.read_snapshot(snapshot);
        nl = md.nlist.cell();

        ############ BONDS ############
        #rigid
        nucl0 = getNucleotide(oligos_list[0][0][0])
        rigid = md.constrain.rigid();
        rigid.set_param('backbone', \
                        types=['sidechain','aux'], \
                        positions = [0.9*nucl0.vectors[0], 0.4*nucl0.vectors[1]]); #magic numbers. Check !!!
        rigid.create_bodies()

        #harmonic
        harmonic1 = md.bond.harmonic()
        harmonic1.bond_coeff.set('backbone', k=5.0, r0=0.75)
        harmonic1.bond_coeff.set('aux_sidechain', k=00.0, r0=0.1) #needed so sidechains in a chain dont interact

        #dihedrals
        def harmonicAngle(theta, kappa, theta0):
           V = 0.5 * kappa * (theta-theta0)**2;
           F = -kappa*(theta-theta0);
           return (V, F)

        index_1st_nucl_in_strand = 0
        for chain in range(len(oligos_list)):
            for strand in range(len(oligos_list[chain])):
                for nucl in range(len(oligos_list[chain][strand]) - 2):

                    bckbone_1 = index_1st_nucl_in_strand + nucl
                    sidechain_1 = total_num_nucl + 2*index_1st_nucl_in_strand + 2*nucl
                    aux_1 = sidechain_1 + 1

                    bckbone_2 = bckbone_1 + 1
                    sidechain_2 = sidechain_1 + 2
                    aux_2 = aux_1 + 2

                    system.dihedrals.add('dihedral1',  sidechain_1, bckbone_1, bckbone_2, sidechain_2)
                    system.dihedrals.add('dihedral21', sidechain_1, bckbone_1, aux_1, bckbone_2)
                    system.dihedrals.add('dihedral22', sidechain_2, bckbone_2, aux_2, bckbone_1)
                    system.dihedrals.add('dihedral31', aux_1, bckbone_1, sidechain_1, bckbone_2)
                    system.dihedrals.add('dihedral32', aux_2, bckbone_2, sidechain_2, bckbone_1)

                    system.bonds.add('aux_sidechain', sidechain_1, sidechain_2)

                index_1st_nucl_in_strand += len(oligos_list[chain][strand])

        dtable = md.dihedral.table(width=1000)
        dtable.dihedral_coeff.set('dihedral1',  func=harmonicAngle, coeff=dict(kappa=50, theta0=-0.28))
        dtable.dihedral_coeff.set('dihedral21', func=harmonicAngle, coeff=dict(kappa=50, theta0=+1.30))
        dtable.dihedral_coeff.set('dihedral22', func=harmonicAngle, coeff=dict(kappa=50, theta0=-1.30))
        dtable.dihedral_coeff.set('dihedral31', func=harmonicAngle, coeff=dict(kappa=50, theta0=-1.57))
        dtable.dihedral_coeff.set('dihedral32', func=harmonicAngle, coeff=dict(kappa=50, theta0=+1.57))

        # fix diameters for vizualization
        for i in range(0, total_num_nucl):
            system.particles[i].diameter = 0.8
        for i in range(total_num_nucl, len(system.particles), 2):
            system.particles[i].diameter = 0.3
            system.particles[i + 1].diameter = 0.1

        ######### INTERACTIONS ############
        # LJ interactions
        wca = md.pair.lj(r_cut=2.0**(1/6), nlist=nl)
        wca.set_params(mode='shift')
        wca.pair_coeff.set('backbone', 'backbone',   epsilon=1.0, sigma=0.750, r_cut=0.750*2**(1/6))
        wca.pair_coeff.set('backbone', 'sidechain',  epsilon=1.0, sigma=0.375, r_cut=0.375*2**(1/6))
        wca.pair_coeff.set('sidechain', 'sidechain', epsilon=1.0, sigma=0.100, r_cut=0.4)
        wca.pair_coeff.set('backbone', 'aux',        epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
        wca.pair_coeff.set('aux', 'sidechain',       epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))
        wca.pair_coeff.set('aux', 'aux',             epsilon=0.0, sigma=1.000, r_cut=1.000*2**(1/6))

        ########## INTEGRATION ############
        md.integrate.mode_standard(dt=0.003, aniso=True);
        rigid = group.rigid_center();
        md.integrate.langevin(group=rigid, kT=0.1, seed=42);
        ########## DUMP & RUN ############
        dump.gsd("output/tripod_MD.gsd",
                       period=1000,
                       group=group.all(),
                       static=[],
                       overwrite=True);
        run(1e6);



#
#
