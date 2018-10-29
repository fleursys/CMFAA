from ComputeTimeStep import computeTimeStep
from InputVariables import nx, ngc, iniCond, bc
from FixedVariables import gamma, nw, rho, v, p, rhov, e
import numpy as np

# returns an array with zeros with nx cells plus ngc ghostcells on each end
def MakeMesh():
    x = [0]*(nx+2*ngc)
    return x

# defines initial conditions for 'shock', the shock tube problem and 'acoustic', the linear acoustic wave problem.
# input is an array x with dimensions nx + 2*ngc, output is the maximum run time, dx and an array x with similar dimensions

def initialization(x):
    if iniCond == 'shock':
        xmin = -0.5
        xmax = 0.5
        tmax = 0.2
        dx = (xmax-xmin)/nx
        for i in range(ngc, ngc + nx//2 + 1):
            prim = [8, 0, 8/gamma]
            x[i] = [rho(prim)]
            x[i].extend(rhov(prim))
            x[i].append(e(prim))
        for i in range(ngc + nx//2 + 1, nx + ngc):
            prim = [1, 0, 1]
            x[i] = [rho(prim)]
            x[i].extend(rhov(prim))
            x[i].append(e(prim))
        return tmax, dx, x
    elif iniCond == 'acoustic':
        xmin = 0
        xmax = 1
        tmax = 3
        dx = (xmax - xmin) / nx
        for i in range(ngc, ngc + nx):
            AcouP = 0.1 + 0.0001*np.exp(-((xmin + (xmax/(nx-1))*(i-2) - 0.5)*(xmin + (xmax/(nx-1))*(i-2) - 0.5))/0.01)
            prim = [AcouP*gamma, 0, AcouP]
            x[i] = [rho(prim)]
            x[i].extend(rhov(prim))
            x[i].append(e(prim))
        return tmax, dx, x
    else:
        print('no valid initial condition')
    return

# defines values for ghost cells for an array x with dimension nx + 2*ngc and returns an array with equal dimensions

def getbc(x):
    if bc == 'fixed':
        for i in range(ngc):
            prim = [8, 0, 8 / gamma]
            x[i] = [rho(prim)]
            x[i].extend(rhov(prim))
            x[i].append(e(prim))
            prim2 = [1, 0, 1]
            x[-(i+1)] = [rho(prim2)]
            x[-(i+1)].extend(rhov(prim2))
            x[-(i+1)].append(e(prim2))
        return x
    elif bc == 'periodic':
        for i in range(ngc):
            x[i] = x[ngc]
            x[-(i+1)] = x[-ngc-1]
        return x
    else:
        print('no valid boundary condition')

# gives the eigenvalues as a 3 by nx+2*ngc matrix lam
def getEigen(x):
    lam = []
    vlam = []
    clam = []
    for i in range(nx + 2*ngc):
        Cs = np.sqrt((gamma*p(x[i])/rho(x[i])))
        lam.append([v(x[i]) - Cs, v(x[i]), v(x[i]) + Cs])
        vlam.append(v(x[i]))
        clam.append(Cs)
    return lam, vlam, clam

# gives the eigenvector matrix for each set of  3 eigenvalues in lam
def getK(lam):
    K = []
    for i in range(nx + 2*ngc):
        K.append([])
        for j in range(nw):
            K[i].append([1, lam[i][j], (lam[i][j]*lam[i][j])/2])
    return K

def getKinv(vlam, clam):
    Kinv = []
    for i in range(nx + 2*ngc):
        Kinv.append([0,0,0])
        Kinv[i][0] = [vlam[i]/(2*clam[i]), -1/2*clam[i]-vlam[i]/(clam[i]*clam[i]), 1/(clam[i]*clam[i])]
        Kinv[i][1] = [1 - vlam[i]*vlam[i]/(clam[i]*clam[i]), 2*vlam[i]/(clam[i]*clam[i]), -2/(clam[i]*clam[i])]
        Kinv[i][2] = [vlam[i]*vlam[i]/(2*clam[i]), -vlam[i]/(clam[i]*clam[i]) + 1/(2*clam[i]), 1/(clam[i]*clam[i])]
    return Kinv

def ConsToDiag(Kinv, x):
    for i in range(nx + 2*ngc):
        w = [0]*nw
        for j in range(nw):
            for k in range(nw):
                w[j] = w[j] + Kinv[i][j][k]*x[i][k]
        x[i] = w
    return x

def DiagtoCons(K, x):
    for i in range(nx + 2*ngc):
        w = [0]*nw
        for j in range(nw):
            for k in range(nw):
                w[j] = w[j] + K[i][j][k]*x[i][k]
        x[i] = w
    return x

def getFlux(lam, x, dx, dt):
    f = [0]*(nx + 2*ngc)
    for l in range(ngc):
        f[l] = [0]*nw
        f[-(l+1)] = [0]*nw
    for i in range(nx):
        f[i+1] = [0]*nw
        for j in range(nw):
            if lam[i+1][j] > 0:
                f[i+1][j] = x[i+1][j] - dt*lam[i+j][j]*((x[i+1][j]-x[i][j])/dx)
            else:
                f[i + 1][j] = x[i + 1][j] - dt * lam[i + j][j] * ((x[i + 2][j] - x[i+1][j]) / dx)
    return f

def IntegrateTime(x,dx):
    lam, vlam, clam = getEigen(x)
    dt = computeTimeStep(lam,dx)
    K = getK(lam)
    Kinv = getKinv(vlam,clam)
    x = ConsToDiag(Kinv, x)
    x = getFlux(lam,x,dx,dt)
    x = DiagtoCons(K, x)
    x = getbc(x)
    return dt, x