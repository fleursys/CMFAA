import numpy as np
from InputVariables import nx, ngc
from FixedVariables import nw


def cmax(lam, nx):
    lmax = 0
    for i in range(nx + 2*ngc):
        for j in range(nw):
            test = np.abs(lam[i][j])
            if test > lmax:
                lmax = test
            else:
                lmax = lmax
    return lmax

def computeTimeStep(lam, dx, nx):
    lamda = cmax(lam, nx)
    return dx/(lamda)

