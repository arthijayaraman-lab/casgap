import numpy as np
from casgap import casgap
params = {
    "seed": 11351,                    # Random Seed
    "boxlength": 100,                 # Length of the box
    "mean_R": 1,                      # Mean of volumetric radii
    "sd_R": 0.001,                    # Standard Deviation of volumetric radii
    "mean_gamma": 2,                  # Mean of aspect ratio
    "sd_gamma": 0.001,                # Standard Deviation of aspect ratio
    "W": [-np.sqrt(6 + 3 * np.sqrt(3)), np.sqrt(2 + np.sqrt(3)), 0] / np.sqrt(8 + 4 * np.sqrt(3)),  # Orientation axis
    "omega": -75,                     # Orientation angle parameter
    "kappa": 1,                       # Orientation anisotropy parameter
    "volfrac": 0.1                    # Target volume fraction
}

restart = {
    "flag": 0,
    "path": './restart/',
    "file": 'session.mat',
    "freq": 1000
}

output = {
    "file": 'sample.dump',
    "path": './output/',
    "freq": 1000
}
casgap(params,restart,output)
