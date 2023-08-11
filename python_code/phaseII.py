from check_overlap import check_overlap
from discretize_ellipsoid import discretize_ellipsoid


def phaseII(particlelist, params):
    if particlelist['Nprime']:
        # Check and adjust coords so that the new ellipsoid does not overlap
        # with any previous ellipsoids.
        successflag, coords = check_overlap(particlelist, params)
        if not successflag:
            return successflag, particlelist

        # Increment population
        population = particlelist['Nprime'] + 1
        particlelist['xyz'][population-1, :] = coords
    else:
        # It is the first particle, there can't be any overlap!
        successflag = 1
        population = 1
        coords = particlelist['xyz'][population-1, :]

    # Append polyhedra
    axislengths = particlelist['ac'][population-1, :]
    quats = particlelist['quat'][population-1, :]
    polyhedron = discretize_ellipsoid(axislengths, quats, coords)
    particlelist['polyhedra'][population-1] = polyhedron
    particlelist['Nprime'] = population
    return successflag, particlelist


