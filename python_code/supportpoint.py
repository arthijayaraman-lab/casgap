import numpy as np
def supportpoint(polyhedron1, polyhedron2, direction):
    furthestpoint1 = furthestpoint(polyhedron1['vertices'], direction)
    furthestpoint2 = furthestpoint(polyhedron2['vertices'], -direction)
    point = furthestpoint1 - furthestpoint2
    return point, furthestpoint1, furthestpoint2


def furthestpoint(vertices, direction):
    maxpoint = np.zeros(3)
    maxprojection = -np.inf
    for i in range(vertices.shape[0]):
        projection = np.dot(vertices[i], direction)
        if maxprojection < projection:
            maxprojection = projection
            maxpoint = vertices[i]
    return maxpoint
