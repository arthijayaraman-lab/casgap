import numpy as np
from scipy.special import erfcinv

def phaseI(params):
    #np.random.seed()  # Reset the random seed

    target_volfrac = params['volfrac']
    box_vol = params['boxlength'] ** 3
    vol_prefactor = 4 / 3 * np.pi
    mean_vol = vol_prefactor * params['mean_R'] ** 3

    # Total number of particles
    N = round(2 * target_volfrac / mean_vol * box_vol)  # Factor of 2 is for good measure.

    # Parameters for lognormal distribution
    R_logmu = np.log(params['mean_R'] ** 2 / np.sqrt(params['mean_R'] ** 2 + params['sd_R'] ** 2))
    gamma_logmu = np.log(params['mean_gamma'] ** 2 / np.sqrt(params['mean_gamma'] ** 2 + params['sd_gamma'] ** 2))
    R_logsigma = np.sqrt(np.log(1 + params['sd_R'] ** 2 / params['mean_R'] ** 2))
    gamma_logsigma = np.sqrt(np.log(1 + params['sd_gamma'] ** 2 / params['mean_gamma'] ** 2))

    # Generate random samples of R
    Ri = lognormrandvar((N, 1), R_logmu, R_logsigma)
    partialvoli = vol_prefactor * Ri ** 3 / box_vol
    actual_volfrac = np.sum(partialvoli)

    while actual_volfrac <= target_volfrac:
        Ri_extra = lognormrandvar((N, 1), R_logmu, R_logsigma)
        Ri = np.concatenate((Ri, Ri_extra))
        voli_extra = vol_prefactor * Ri_extra ** 3 / box_vol
        partialvoli = np.concatenate((partialvoli, voli_extra))
        actual_volfrac = np.sum(partialvoli)

    if actual_volfrac > target_volfrac:
        N = np.argmax(np.cumsum(partialvoli) > target_volfrac)  # Find index where cumulative sum exceeds target_volfrac
        Ri = Ri[:N]

    # Generate random samples of gamma
    gammai = lognormrandvar((N, 1), gamma_logmu, gamma_logsigma)

    # Set spheroid axis lengths
    ai = Ri / gammai ** (1 / 3)
    ci = Ri * gammai ** (2 / 3)
    ac = np.hstack((ai, ci))
    Lambda = generate_quat(params['W'], params['omega'])
    quat = samplequat(N, Lambda, params['kappa'])

    particlelist = {'N': N, 'ac': ac, 'quat': quat}

    print(f"Phase I completed! N is {N}.")
    return particlelist


def lognormrandvar(size, logmu, logsigma):
    temp = -np.sqrt(2) * erfcinv(2 * np.random.rand(*size))  # Std normal random variable
    result = np.exp(logmu + temp * logsigma)
    return result


def generate_quat(W, omega):
    alpha_by_2 = omega / 2 * np.pi / 180  # Convert angle to radians
    norm = W / np.sqrt(np.sum(W** 2))
    q = np.array([np.cos(alpha_by_2),norm[0] * np.sin(alpha_by_2),norm[1] * np.sin(alpha_by_2),norm[2] * np.sin(alpha_by_2)])
    return q


def samplequat(num, pref_q, kappa):
    v = -np.sqrt(2) * erfcinv(2 * np.random.rand(num, 2))
    vsum = np.sqrt(np.sum(v ** 2, axis=1))
    v = v / vsum[:, np.newaxis]  # Ensure shape compatibility

    if kappa:
        randnum = np.random.rand(num)
        w = 1 + 1 / kappa * np.log(randnum + (1 - randnum) * np.exp(-2 * kappa))
    else:
        w = 2 * np.random.rand(num) - 1

    orgvec = np.hstack((np.zeros((num, 2)), np.ones((num, 1))))
    newvec = np.vstack((w, np.sqrt(1 - w ** 2) * v[:, 0], np.sqrt(1 - w ** 2) * v[:, 1])).T

    axisvec = np.cross(orgvec, newvec)
    axisvec /= np.sqrt(np.sum(axisvec ** 2, axis=1, keepdims=True))  # Ensure shape compatibility
    axistheta = np.arccos(np.sum(orgvec * newvec, axis=1))

    cos_half = np.cos(axistheta / 2)[:, np.newaxis]  # Reshape to (num, 1)
    sin_half = np.sin(axistheta / 2)[:, np.newaxis]  # Reshape to (num, 1)
    quats = np.hstack((cos_half, sin_half * axisvec))  # Concatenate along axis 1

    if kappa:
        Mumat = np.zeros((4, 4))
        Mumat[:, 0] = pref_q.T
        Q, R = np.linalg.qr(Mumat)
        if R[0, 0] < 0:
            Q = -Q
        result = np.dot(Q, quats.T).T
    else:
        result = quats

    return result


