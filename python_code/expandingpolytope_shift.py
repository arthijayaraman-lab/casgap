import numpy as np
from supportpoint import supportpoint
def expandingpolytope_shift(polyhedron1, polyhedron2, simplex):
    polytope = {
        'vertices': simplex['vertices'],
        'faces': simplex['faces']
    }
    facenorms = getfacenormals(polytope)
    tolerance = 1e-6
    #counter=0
    while True:
        #counter += 1
        #print(f"I am stuck here. {counter}")
        minInd, mindist = findnearestface(polytope, facenorms)
        nearestnormal = facenorms[minInd]
        support = supportpoint(polyhedron1, polyhedron2, nearestnormal)[0]
        supportdist = np.dot(nearestnormal, support)
        if (supportdist - mindist) > tolerance:
            numnormals = facenorms.shape[0]
            uniqueedges = []
            newfaces = []
            newnormals = []
            for i in range(numnormals):
                face = polytope['faces'][i]
                facevertex1 = polytope['vertices'][face[0]-1]
                if np.dot(facenorms[i], support - facevertex1) > tolerance:
                    faceedges = np.column_stack((face[:-1], face[1:]))
                    faceedges = np.vstack((faceedges, [face[-1], face[0]]))
                    uniqueedges = finduniqueedges(uniqueedges, faceedges)
                else:
                    newfaces.append(polytope['faces'][i])
                    newnormals.append(facenorms[i])
            # add the new support point to the polytope
            polytope['vertices'] = np.vstack((polytope['vertices'], support))
            supportind = polytope['vertices'].shape[0]
            numuniqueedges = uniqueedges.shape[0]
            for i in range(numuniqueedges):
                newaddedface, newaddedfacenorm = createnewpolytopeface(polytope, uniqueedges[i], supportind)
                newfaces.append(newaddedface)
                newnormals.append(newaddedfacenorm)
            polytope['faces'] = newfaces
            facenorms = np.array(newnormals)
        else:
            shiftvector = nearestnormal
            shiftdist = mindist
            break

    return shiftvector, shiftdist

def getfacenormals(polytope):
    faces = polytope['faces']
    numfaces = len(faces)
    facenorms = np.zeros((numfaces, 3))
    tolerance = 1e-6
    for i in range(numfaces):
        face = faces[i]
        point_a = polytope['vertices'][face[0]-1]
        point_b = polytope['vertices'][face[1]-1]
        point_c = polytope['vertices'][face[2]-1]
        vec_ab = point_b - point_a
        vec_bc = point_c - point_b
        facenorm_abc = np.cross(vec_ab, vec_bc)
        facedist = np.dot(point_a, facenorm_abc)
        if facedist < tolerance:
            facenorm_abc = -facenorm_abc
        facenorms[i] = facenorm_abc / np.linalg.norm(facenorm_abc)  # These are unit vectors
               
    return facenorms


def finduniqueedges(uniqueedges, faceedges):
    if len(uniqueedges) == 0:
        uniqueedges = faceedges
    else:
        for i in range(faceedges.shape[0]):
            # check if the reverse of  uniqueedge is in the faceedge. If there is, then delete.
            uniqueedge_vert1 = uniqueedges[:, 0] - faceedges[i, 1]
            uniqueedge_vert2 = uniqueedges[:, 1] - faceedges[i, 0]
            uniqueness_check = np.logical_not(uniqueedge_vert1) & np.logical_not(uniqueedge_vert2)
            if np.any(uniqueness_check):
                uniqueedges = uniqueedges[np.logical_not(uniqueness_check), :]
            else:
                uniqueedges = np.vstack((uniqueedges, faceedges[i, :]))
    return uniqueedges


def createnewpolytopeface(polytope, edge, pointind):
    point_a = polytope['vertices'][edge[0]-1, :]
    point_b = polytope['vertices'][edge[1]-1, :]
    point_c = polytope['vertices'][pointind-1, :]
    tolerance = 1e-6
    vec_ab = point_b - point_a
    vec_bc = point_c - point_b
    facenorm_abc = np.cross(vec_ab, vec_bc)
    facenorm = facenorm_abc / np.linalg.norm(facenorm_abc)  # These are unit vectors
    facedist = np.dot(point_a, facenorm)
    if facedist > tolerance:
        face = [edge[0], edge[1], pointind]
    else:
        face = [edge[1], edge[0], pointind]
        facenorm = -facenorm
    return face, facenorm

def findnearestface(polytope, facenorms):
    #distances = np.dot(polytope['vertices'][x[0]-1 for x in polytope['faces']], facenorms.T)
    facevertex=np.array([polytope['vertices'][x[0]-1] for x in polytope['faces']])
    distances = np.einsum("ij,ij->i",facevertex,facenorms)
    minInd = np.argmin(distances)
    mindist = distances[minInd]
    return minInd, mindist

