import numpy as np
from InputVariables import nx, ngc

def lamMax(lam):
    lmax = 0
    for i in range(nx + 2*ngc):
        test = np.abs(lam[i][2])
        if test > lmax:
            lmax = test
        else:
            lmax = lmax
    return lmax

def computeTimeStep(lam,dx):
    lmax = lamMax(lam)
    return dx/lmax

