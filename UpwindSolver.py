import numpy as np
# variables to be changed
nx = 10
ngc = 2
iniCond = 'acoustic'
# iniCond = 'shock'
# bc = 'fixed'
bc = 'periodic'

#variables that are fixed
gamma = 1.4

#conservative variables in an array w with dimensions 3,1
def rho(w):
    return w[0]

def rhov(w):
    return w[1]

def e(w):
    return w[2]

def v(w):
    return rhov(w)/rho(w)

def p(w):
    return (gamma-1)*(-rhov(w)*rhov(w))/(2*rho(w))+ e(w)


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
        for i in range(ngc, ngc + nx//2):
            x[i] = [8, 0, 8/gamma]
        for i in range(ngc + nx//2 + 1, nx + ngc):
            x[i] = [1, 0, 1]
        return xmin, xmax, x
    elif iniCond == 'acoustic':
        xmin = 0
        xmax = 1
        for i in range(ngc, ngc + nx):
            AcouP = 0.1 + 0.0001*np.exp(-((xmin + (xmax/(nx-1))*(i-2) - 0.5)*(xmin + (xmax/(nx-1))*(i-2) - 0.5))/0.01)
            x[i] = [AcouP*gamma, 0, AcouP]
        return xmin, xmax, x
    else:
        print('no valid initial condition')
    return

# defines values for ghost cells for an array x with dimension nx + 2*ngc and returns an array with equal dimensions

def getbc(x):
    if bc == 'fixed':
        for i in range(ngc):
            x[i] = [8, 0, 8/gamma]
            x[i-1] = [1, 0, 1]
        return x
    elif bc == 'periodic':
        for i in range(ngc):
            x[i] = x[ngc]
            x[-(i+1)] = x[-ngc-1]
        return x
    else:
        print('no valid boundary condition')

# gives the eigenvalues
def getEigen(x):
    lam = []
    for i in range(len(x)):
        Cs = np.sqrt((gamma*p(x[i])/rho(x[i])))
        lam.append([v(x[i]) - Cs, v(x[i]), v(x[i]) + Cs])
    return lam

def getK(lam):

def getKinv(lam):

def IntegrateTime():
    mesh = MakeMesh()
    xmin, xmax, mesh = initialization(mesh)
    mesh = getbc(mesh)
    eigenval = getEigen(mesh)
    return eigenval

print(IntegrateTime())
