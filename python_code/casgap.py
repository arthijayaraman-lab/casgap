import os
import numpy as np
import pickle
from PhaseI import phaseI
from phaseII import phaseII
from write_data import write_data
import time


def casgap(params, restart, output):
    # Initialize random number generator
    np.random.seed(params['seed'])
    currentseed = np.random.default_rng()
    os.makedirs(restart['path'], exist_ok=True)
    os.makedirs(output['path'], exist_ok=True)

    if restart['flag']:
        with open(os.path.join(restart['path'], restart['file']), 'rb') as f:
            particlelist, loopstart, loopend = pickle.load(f)
        np.random.default_rng(currentseed)
    else:
        particlelist = phaseI(params)
        #particlelist['polyhedra'] = {'vertices': [], 'faces': {}, 'center': []}
        particlelist['polyhedra'] = {}
        particlelist['Nprime'] = 0
        particlelist['xyz'] = (np.random.rand(particlelist['N'], 3) - 0.5) * params['boxlength']
        loopstart = 1
        loopend = particlelist['N']
        currentseed = np.random.default_rng()
        with open(os.path.join(restart['path'], restart['file']), 'wb') as f:
            pickle.dump((particlelist, loopstart, loopend), f)

    start_time = time.time()
    for i in range(loopstart, loopend + 1):
        successflag, particlelist = phaseII(particlelist, params)
        if not successflag:
            break

        if i % restart['freq'] == 0:
            currentseed = np.random.default_rng()
            loopstart = i
            with open(os.path.join(restart['path'], restart['file']), 'wb') as f:
                pickle.dump((particlelist, loopstart, loopend), f)
            end_time1 = time.time()
            print(f"Writing Restart File! Current population is {i}. The time elapsed is {end_time1-start_time} seconds.")

        if i % output['freq'] == 0:
            write_data(output, particlelist, params)

    currentseed = np.random.default_rng()
    with open(os.path.join(restart['path'], restart['file']), 'wb') as f:
        pickle.dump((particlelist, loopstart, loopend), f)
    write_data(output, particlelist, params)

    partialvoli = 4 / 3 * np.pi * particlelist['ac'][:particlelist['Nprime'], 0] ** 2 * particlelist['ac'][
                                                                                        :particlelist['Nprime'], 1]
    actual_volfrac = np.sum(partialvoli) / params['boxlength'] ** 3
    end_time2 = time.time()
    if not successflag:
        print(
            f"Not all ellipsoids could be added. Final population is {particlelist['Nprime']}. Final volume fraction is {actual_volfrac}. The time elapsed is {end_time2 - start_time} seconds.")
    else:
        print(
            f"All ellipsoids were successfully added! Final population is {particlelist['Nprime']}. Final volume fraction is {actual_volfrac}. The time elapsed is {end_time2 - start_time} seconds.")
