import numpy as np
from timeit import default_timer as timer

start = timer()

# variables to be changed
nx = 10
ngc = 1
#iniCond = 'acoustic'
iniCond = 'shock'
bc = 'fixed'
#bc = 'periodic'

# variables that are fixed
gamma = 1.4
nw = 3

# conservative variables in an array w with dimensions 3,1
def rho(w):
    return w[0]

def rhov(w):
    return w[1]

def e(w):
    return w[-1]

def v(w):
    return rhov(w)/rho(w)

def p(w):
    return (gamma-1)*(-rhov(w)*rhov(w))/(2*rho(w)) + e(w)


# returns an array with zeros with nx cells plus ngc ghostcells on each end
def MakeMesh():
    x = [0]*(nx+2*ngc)
    return x

# defines initial conditions for 'shock', the shock tube problem and 'acoustic, the linear acoustic wave problem.
# input is an array x with dimensions nx + 2*ngc, output is gamma, boundaries of the x domain, and an array x with
# the same dimensions as the input

def initialization(x):
    if iniCond == 'shock':
        xmin = -0.5
        xmax = 0.5
        tmax = 0.2
        dx = (xmax-xmin)/nx
        for i in range(ngc, ngc + nx//2 + 1):
            x[i] = [8, 0, 8/gamma]
        for i in range(ngc + nx//2 + 1, nx + ngc):
            x[i] = [1, 0, 1]
        return tmax, dx, x
    elif iniCond == 'acoustic':
        xmin = 0
        xmax = 1
        tmax = 3
        dx = (xmax - xmin) / nx
        for i in range(ngc, ngc + nx):
            AcouP = 0.1 + 0.0001*np.exp(-((xmin + (xmax/(nx-1))*(i-2) - 0.5)*(xmin + (xmax/(nx-1))*(i-2) - 0.5))/0.01)
            x[i] = [AcouP*gamma, 0, AcouP]
        return tmax, dx, x
    else:
        print('no valid initial condition')
    return

# defines values for ghost cells for an array x with dimension nx + 2*ngc and returns an array with equal dimensions

def getbc(x):
    if bc == 'fixed':
        for i in range(ngc):
            x[i] = [8, 0, 8/gamma]
            x[-(i+1)] = [1, 0, 1]
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

def UpwindSolver():
    x = MakeMesh()
    tmax, dx, x = initialization(x)
    x = getbc(x)
    t = 0
    while t < tmax:
        dt, x = IntegrateTime(x,dx)
        t += dt
    return

UpwindSolver()

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ")
print(chrono(start,end))
