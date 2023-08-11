import numpy as np
from supportpoint import supportpoint

def gjk_simplex(polyhedron1, polyhedron2):
    direction = polyhedron1['center'] - polyhedron2['center']
    support = supportpoint(polyhedron1, polyhedron2, direction)[0]
    simplexvertices = support
    direction = -np.array(support)

    while True:
        support = supportpoint(polyhedron1, polyhedron2, direction)[0]

        if np.all(np.logical_not(direction)) or np.dot(support, direction) <= 0:
            intersection_flag = 0
            simplex = {'vertices': []}
            break

        simplexvertices = np.vstack((support, simplexvertices))
        flag, simplexvertices, direction = nextsimplex(simplexvertices)

        if flag:
            if simplexvertices.shape[0] < 4:
                intersection_flag = 0
                simplex = {'vertices': []}
            else:
                intersection_flag = 1
                simplexfaces = gettetrahedronfaces(simplexvertices, [0, 0, 0])
                simplex = {'vertices': simplexvertices, 'faces': simplexfaces}
            break

    return intersection_flag, simplex


def nextsimplex(vertices):
    num_vertices = len(vertices)
    if num_vertices == 2:
        flag, vertices, direction = linesimplex(vertices)
    elif num_vertices == 3:
        flag, vertices, direction = trianglesimplex(vertices)
    elif num_vertices == 4:
        flag, vertices, direction = tetrahedronsimplex(vertices)
    else:
        flag = False
        vertices = []
        direction = np.array([])

    return flag, vertices, direction

def linesimplex(vertices):
    direction = np.zeros(3)
    point_a = vertices[0]
    point_b = vertices[1]
    vec_ab = point_b - point_a
    vec_ao = -point_a
    norm_abo = np.cross(vec_ab, vec_ao)

    if np.linalg.norm(norm_abo):  # check for collinearity
        direction = np.cross(norm_abo, vec_ab)
        flag = False
    else:
        flag = True

    return flag, vertices, direction

def trianglesimplex(vertices):
    direction = np.zeros(3)
    point_a = vertices[0]
    point_b = vertices[1]
    point_c = vertices[2]

    vec_ab = point_b - point_a
    vec_ac = point_c - point_a
    vec_ao = -point_a

    facenorm_abc = np.cross(vec_ab, vec_ac)

    if np.linalg.norm(facenorm_abc):  # check if origin is coplanar
        if np.dot(facenorm_abc, vec_ao) > 0:
            direction = facenorm_abc
            flag = False
        else:
            direction = -facenorm_abc
            flag = False
    else:
        flag = True

    return flag, vertices, direction

def tetrahedronsimplex(vertices):
    direction = np.zeros(3)
    point_a = vertices[0]
    point_b = vertices[1]
    point_c = vertices[2]
    point_d = vertices[3]

    vec_ab = point_b - point_a
    vec_ac = point_c - point_a
    vec_ad = point_d - point_a
    vec_ao = -point_a

    facenorm_abc = np.cross(vec_ab, vec_ac)
    if np.dot(facenorm_abc, vec_ad) > 0:
        facenorm_abc = -facenorm_abc

    facenorm_acd = np.cross(vec_ac, vec_ad)
    if np.dot(facenorm_acd, vec_ab) > 0:
        facenorm_acd = -facenorm_acd

    facenorm_abd = np.cross(vec_ab, vec_ad)
    if np.dot(facenorm_abd, vec_ac) > 0:
        facenorm_abd = -facenorm_abd

    outsideabc = np.dot(facenorm_abc, vec_ao) > 0
    outsideacd = np.dot(facenorm_acd, vec_ao) > 0
    outsideabd = np.dot(facenorm_abd, vec_ao) > 0

    if (outsideabc and outsideacd) or (outsideacd and outsideabd) or (outsideabd and outsideabc):
        flag = False
    elif outsideabc:
        vertices = np.array([point_a, point_b, point_c])
        direction = facenorm_abc
        flag = False
    elif outsideacd:
        vertices = np.array([point_a, point_c, point_d])
        direction = facenorm_acd
        flag = False
    elif outsideabd:
        vertices = np.array([point_a, point_b, point_d])
        direction = facenorm_abd
        flag = False
    else:
        flag = True

    return flag, vertices, direction


def gettetrahedronfaces(vertices, interiorpoint):
    if vertices.shape[0] != 4:
        raise ValueError('Input should be a set of tetrahedral vertices')

    faces = [[1, 2, 3], [1, 4, 2], [1, 3, 4], [2, 4, 3]]

    for i in range(4):
        face = faces[i]
        point_a = vertices[face[0] - 1]
        point_b = vertices[face[1] - 1]
        point_c = vertices[face[2] - 1]

        vec_ab = point_b - point_a
        vec_bc = point_c - point_b
        facenorm_abc = np.cross(vec_ab, vec_bc)

        if np.dot(facenorm_abc, point_b - interiorpoint) < 0:
            faces[i] = [face[0], face[2], face[1]]

    return faces





