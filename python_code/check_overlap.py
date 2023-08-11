import numpy as np
from discretize_ellipsoid import discretize_ellipsoid
from gjk_simplex import gjk_simplex
from expandingpolytope_shift import expandingpolytope_shift

def check_overlap(particlelist, params):
    successflag = 0

    population = particlelist['Nprime']
    newind = population + 1
    ac = particlelist['ac'][newind - 1]
    q = particlelist['quat'][newind - 1]
    boundingradius = np.max(ac)
    newpolyhedra_unshifted = discretize_ellipsoid(ac, q, [0, 0, 0])

    maxattempts = 10000
    perturbfreq = 100  # resample new coords at this rate
    maxradii = np.max(particlelist['ac'][:population], axis=1)
    pos = particlelist['xyz'][:population]

    EPAclearance = 0.1
    for numattempts in range(1, maxattempts + 1):
        if numattempts % perturbfreq == 1:
            if numattempts < perturbfreq:
                coords = (np.random.rand(3) - 0.5) * params['boxlength']
            else:
                randvars = np.random.rand(3)
                zbins = np.linspace(-0.5 * params['boxlength'], 0.5 * params['boxlength'], 11)
                znum, _ = np.histogram(pos[:, 2], bins=zbins)
                maxznum = np.max(znum)
                minznum = np.min(znum)
                if maxznum == minznum:
                    zcoord = (randvars[0] - 0.5) * params['boxlength']
                else:
                    newznum = (maxznum - znum) / (maxznum - minznum)
                    cumznum = np.concatenate(([0],np.cumsum(newznum))) + np.arange(11) * 1e-10
                    cumznum /= cumznum[-1]
                    zcoord = np.interp(randvars[2], cumznum, zbins)
                zbinwidth = zbins[1] - zbins[0]
                zfilter = (pos[:, 2] > zcoord - zbinwidth / 2) & (pos[:, 2] <= zcoord + zbinwidth / 2)

                ybins = np.linspace(-0.5 * params['boxlength'], 0.5 * params['boxlength'], 11)
                ynum, _ = np.histogram(pos[zfilter, 1], bins=ybins)
                maxynum = np.max(ynum)
                minynum = np.min(ynum)
                if not np.sum(zfilter) or maxynum == minynum:
                    ycoord = (randvars[1] - 0.5) * params['boxlength']
                else:
                    newynum = (maxynum - ynum) / (maxynum - minynum)
                    cumynum = np.concatenate(([0],np.cumsum(newynum))) + np.arange(11) * 1e-10
                    cumynum /= cumynum[-1]
                    ycoord = np.interp(randvars[1], cumynum, ybins)
                ybinwidth = ybins[1] - ybins[0]
                yzfilter = zfilter & (pos[:, 1] > ycoord - ybinwidth / 2) & (pos[:, 1] <= ycoord + ybinwidth / 2)

                xbins = np.linspace(-0.5 * params['boxlength'], 0.5 * params['boxlength'], 11)
                xnum, _ = np.histogram(pos[yzfilter, 0], bins=xbins)
                maxxnum = np.max(xnum)
                minxnum = np.min(xnum)
                if not np.sum(yzfilter) or maxxnum == minxnum:
                    xcoord = (randvars[0] - 0.5) * params['boxlength']
                else:
                    newxnum = (maxxnum - xnum) / (maxxnum - minxnum)
                    cumxnum = np.concatenate(([0],np.cumsum(newxnum))) + np.arange(11) * 1e-10
                    cumxnum /= cumxnum[-1]
                    xcoord = np.interp(randvars[0], cumxnum, xbins)
                coords = np.array([xcoord, ycoord, zcoord])

        newpolyhedra = newpolyhedra_unshifted.copy()
        newpolyhedra['vertices'] = newpolyhedra['vertices'] + np.ones((newpolyhedra['vertices'].shape[0], 1)) * coords
        dist = np.sqrt(np.sum((np.ones((population, 1)) * coords - pos) ** 2, axis=1))
        potentialoverlap = (maxradii + boundingradius - dist > 0)

        if not np.any(potentialoverlap):
            successflag = 1
            break

        indices = np.where(potentialoverlap)[0]
        numintersections = 0
        cummulativeshift = np.array([0, 0, 0])
        for i in indices:
            #try:
            intersectionflag, intersectionsimplex = gjk_simplex(newpolyhedra, particlelist['polyhedra'][i])
            #except KeyError:
            #    continue

            if intersectionflag:
                numintersections += 1
                shiftvector, shiftdist = expandingpolytope_shift(newpolyhedra, particlelist['polyhedra'][i],
                                                                 intersectionsimplex)
                cummulativeshift = -shiftvector * (shiftdist + EPAclearance)

        if np.any(cummulativeshift):
            if not numattempts % 100:
                print(
                    f"Warning: Previous attempts failed to converge after {numattempts} attempts. Current ellipsoid has {numintersections} intersections. The population is {population}.")
            coords = coords + cummulativeshift
            coords = np.mod(coords + params['boxlength'] / 2, params['boxlength']) - params['boxlength'] / 2
        else:
            successflag = 1
            break

    return successflag, coords

