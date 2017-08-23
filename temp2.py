def distanceBetweenVhs(vh1, index1, vh2, index2):
    pos1 = findCoordinates(vh1, index1)[1]
    pos2 = findCoordinates(vh2, index2)[1]
    distance = np.linalg.norm(pos1 - pos2)
    return(distance)

def connection3p(strand):
    if strand.connection3p() != None:
            vh1 = strand.idNum()
            index1 = strand.connection3p().idx5Prime()
            vh2 = strand.connection3p().idNum()
            index2 = strand.connection3p().connection5p().idx3Prime()
            distance = distanceBetweenVhs(vh1, index1, vh2, index2)
            if distance < 10.0:
                return(vh2)

def connection5p(strand):
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
    This function separates the origami 'part' into bodies
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











#t
