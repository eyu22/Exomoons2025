import pymultinest
import numpy as np

def prior(cube):
    cube[0] = cube[0] * 10
    return cube

def loglike(cube):
    return -0.5 * ((cube[0] - 5.0) ** 2)

pymultinest.run(loglike, prior, 1, outputfiles_basename='test_', resume=False, verbose=True)
