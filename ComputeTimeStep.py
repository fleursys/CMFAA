import numpy as np
from InputVariables import nx, ngc
from FixedVariables import nw

def lamMax(lam):
    lmax = 0.000000000001
    for i in range(nx + 2*ngc):
        for j in range(nw):
            test = np.abs(lam[i][j])
            if test > lmax:
                lmax = test
            else:
                lmax = lmax
    return lmax

def computeTimeStep(lam,dx):
    lmax = lamMax(lam)
    return dx/lmax

