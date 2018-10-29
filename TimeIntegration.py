from ComputeTimeStep import computeTimeStep
from InputVariables import nx, ngc, iniCond, bc
from FixedVariables import gamma, nw, rho, v, p, rhov, e
import numpy as np

# returns an array with zeros with nx cells plus ngc ghostcells on each end
def Grid():
    x = [0]*(nx+2*ngc)
    return x

# defines initial conditions for 'shock', the shock tube problem and 'acoustic', the linear acoustic wave problem.
# input is an array x with dimensions nx + 2*ngc, output is the maximum run time, dx and an array x with similar dimensions

def initialization():
    x = Grid()
    if iniCond == 'shock':
        xmin = -0.5
        xmax = 0.5
        tmax = 0.2
        dx = (xmax-xmin)/(nx-1)
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
        return tmax, xmax, xmin, dx, x
    elif iniCond == 'acoustic':
        xmin = 0
        xmax = 1
        tmax = 3
        dx = (xmax - xmin) / (nx-1)
        for i in range(ngc, ngc + nx):
            AcouP = 0.1 + 0.0001*np.exp(-((xmin + (xmax/(nx-1))*(i-2) - 0.5)*(xmin + (xmax/(nx-1))*(i-2) - 0.5))/0.01)
            prim = [AcouP*gamma, 0, AcouP]
            x[i] = [rho(prim)]
            x[i].extend(rhov(prim))
            x[i].append(e(prim))
        return tmax, xmax, xmin, dx, x
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

# gives the eigenvalues as a 3 by nx+2*ngc matrix lam, it also gives the velocity (vlam) and the speed of sound (clam)
def getEigen(x):
    lam = []
    vlam = []
    clam = []
    for i in range(nx + 2*ngc):
        Cs = np.sqrt((gamma*p(x[i])/rho(x[i])))
        lam.append([v(x[i])[0] - Cs, v(x[i])[0], v(x[i])[0] + Cs])
        vlam.append(v(x[i])[0])
        clam.append(Cs)
    return lam, vlam, clam

# gives the eigenvector matrix for each point in the mesh based on the the eigenvalue array lam
def getK(lam):
    K = []
    for i in range(nx + 2*ngc):
        K.append([])
        for j in range(nw):
            K[i].append([1, lam[i][j], (lam[i][j]*lam[i][j])/2])
    return K

# gives the inverse eigenvector matrix for each point in the mesh for the velocity and the speed of sound
def getKinv(vlam, clam):
    Kinv = Grid()
    for i in range(nx + 2*ngc):
        Kinv[i] = [0, 0, 0]
        Kinv[i][0] = [vlam[i]/(2*clam[i]) + vlam[i]*vlam[i]/(2*clam[i]*clam[i]), -1/(2*clam[i]) - vlam[i]/(clam[i]*clam[i]), 1/(clam[i]*clam[i])]
        Kinv[i][1] = [1 - vlam[i]*vlam[i]/(clam[i]*clam[i]), 2*vlam[i]/(clam[i]*clam[i]), -2/(clam[i]*clam[i])]
        Kinv[i][2] = [vlam[i]*vlam[i]/(2*clam[i]) - vlam[i]/(2*clam[i]), -vlam[i]/(clam[i]*clam[i]) + 1/(2*clam[i]), 1/(clam[i]*clam[i])]
    return Kinv

# gives the multiplication of K-1 with the conservative variables in each point, resulting in the diagonal variables
def ConsToDiag(Kinv, x):
    f = Grid()
    for i in range(nx + 2*ngc):
        f[i] = []
        for j in range(nw):
            f[i].append(Kinv[i][j][0]*x[i][0] + Kinv[i][j][1]*x[i][1] + Kinv[i][j][2]*x[i][2])
    return f

# gives the multiplication of K with the diagonal variables in each point, resulting in the conservation variables
def DiagtoCons(K, x):
    f = Grid()
    for i in range(nx + 2*ngc):
        f[i] = []
        for j in range(nw):
            f[i].append(K[i][0][j]*x[i][0] + K[i][1][j]*x[i][1] + K[i][2][j]*x[i][2])
    return f

# computes the diagonal variables in the next time step (after dt). The ghostcells are given as zero.
def getFlux(lam, x, dx, dt):
    f = Grid()
    for l in range(ngc):
        f[l] = [0]*nw
        f[-(l+1)] = [0]*nw
    for i in range(nx):
        f[i+1] = [0]*nw
        for j in range(nw):
            if lam[i+1][j] > 0:
                f[i + 1][j] = x[i + 1][j] - dt/dx * lam[i + 1][j] * (x[i + 1][j] - x[i][j])
            else:
                f[i + 1][j] = x[i + 1][j] - dt/dx * lam[i + 1][j] * (x[i + 2][j] - x[i + 1][j])
    return f

# computes the conservative variables in the next time step (after dt).
def IntegrateTime(x, dx):
    x = getbc(x)
    lam, vlam, clam = getEigen(x)
    dt = computeTimeStep(lam, dx)
    K = getK(lam)
    Kinv = getKinv(vlam, clam)
    diag = ConsToDiag(Kinv, x)
    flux = getFlux(lam, diag, dx, dt)
    cons = DiagtoCons(K, flux)
    return dt, cons