import numpy as np
def discretize_ellipsoid(axislengths, q, origin):
    # Find the axial vector from the quaternion
    axvec = np.array([
        [q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2, 2*(q[1]*q[2] - q[0]*q[3]), 2*(q[1]*q[3] + q[0]*q[2])],
        [2*(q[1]*q[2] + q[0]*q[3]), q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2, 2*(q[2]*q[3] - q[0]*q[1])],
        [2*(q[1]*q[3] - q[0]*q[2]), 2*(q[2]*q[3] + q[0]*q[1]), q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2]
    ])
    mag_axvec = np.sqrt(axvec[0]**2 + axvec[1]**2 + axvec[2]**2)
    axvec = axvec / mag_axvec

    # Convert Ellipsoid to Polyhedron
    f = np.sqrt(2) - 1
    a, b, c = axislengths[0], axislengths[0], axislengths[1]
    X0 = a * np.array([1, f, f])
    Y0 = b * np.array([f, 1, f])
    Z0 = c * np.array([f, f, 1])
    allvertices = np.concatenate([
        np.column_stack((X0, Y0, Z0)),  # first octant: x y z
        np.column_stack((-X0, Y0, Z0)),  # second octant: -x y z
        np.column_stack((-X0, -Y0, Z0)),  # third octant: -x -y z
        np.column_stack((X0, -Y0, Z0)),  # fourth octant: x -y z
        np.column_stack((X0, Y0, -Z0)),  # fifth octant: x y -z
        np.column_stack((-X0, Y0, -Z0)),  # sixth octant: -x y -z
        np.column_stack((-X0, -Y0, -Z0)),  # seventh octant: -x -y -z
        np.column_stack((X0, -Y0, -Z0))  # eighth octant: x -y -z
    ])
    numvertices = allvertices.shape[0]
    allvertices = np.dot(allvertices, axvec.T) + np.ones((numvertices, 1)) * origin
    allfaces = [
        [1, 10, 22, 13],  # +x rect
        [4, 16, 19, 7],  # -x rect
        [2, 14, 17, 5],  # +y rect
        [8, 20, 23, 11],  # -y rect
        [3, 6, 9, 12],  # +z rect
        [15, 24, 21, 18],  # -z rect
        [1, 3, 12, 10],  # +x+z rect
        [4, 7, 9, 6],  # -x+z rect
        [16, 18, 21, 19],  # -x-z rect
        [13, 22, 24, 15],  # +x-z rect
        [2, 5, 6, 3],  # +y+z rect
        [8, 11, 12, 9],  # -y+z rect
        [20, 21, 24, 23],  # -y-z rect
        [14, 15, 18, 17],  # +y-z rect
        [1, 13, 14, 2],  # +x+y rect
        [4, 5, 17, 16],  # -x+y rect
        [7, 19, 20, 8],  # -x-y rect
        [10, 11, 23, 22],  # +x-y rect
        [1, 2, 3],  # +x+y+z tri
        [4, 6, 5],  # -x+y+z tri
        [7, 8, 9],  # -x-y+z tri
        [10, 12, 11],  # +x-y+z tri
        [13, 15, 14],  # +x+y-z tri
        [16, 17, 18],  # -x+y-z tri
        [19, 21, 20],  # -x-y-z tri
        [22, 23, 24]  # +x-y-z tri
    ]
    polyhedra = {
        'vertices': allvertices,
        'faces': allfaces,
        'center': origin
    }

    return polyhedra
