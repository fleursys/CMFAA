import numpy as np
# variables to be changed
nx = 10
ngc = 2
iniCond = 'acoustic'
# iniCond = 'shock'
# bc = 'fixed'
bc = 'periodic'

# returns an array with zeros with nx cells plus ngc ghostcells on each end
def MakeMesh():
    mesh = [0]*(nx+2*ngc)
    return mesh

# defines initial conditions for 'shock', the shock tube problem and 'acoustic, the linear acoustic wave problem.
def initialization(mesh):
    if iniCond == 'shock':
        gamma = 1.4
        xmin = -0.5
        xmax = 0.5
        for i in range(ngc, ngc + nx//2):
            mesh[i] = [8, 0, 8/gamma]
        for i in range(ngc + nx//2 + 1, nx + ngc):
            mesh[i] = [1, 0, 1]
        return gamma, xmin, xmax, mesh
    elif iniCond == 'acoustic':
        gamma = 1.4
        xmin = 0
        xmax = 1
        for i in range(ngc, ngc + nx):
            AcouP = 0.1 + 0.0001*np.exp(-((xmin + (xmax/(nx-1))*(i-2) - 0.5)*(xmin + (xmax/(nx-1))*(i-2) - 0.5))/0.01)
            mesh[i] = [AcouP*gamma, 0, AcouP]
        return gamma, xmin, xmax, mesh
    else:
        print('no valid initial condition')
    return

# defines values for ghost cells
def getbc(mesh):
    if bc == 'fixed':
        for i in range(ngc):
            mesh[i] = [8, 0, 8/gamma]
            mesh[i-1] = [1, 0, 1]
    elif bc == 'periodic':
        for i in range(ngc):
            mesh[i] = mesh[ngc]
            mesh[-(i+1)] = mesh[-ngc-1]
        return mesh
    else:
        print('no valid boundary condition')

def UpwindSolver():
    mesh = MakeMesh()
    gamma, xmin, xmax, mesh = initialization(mesh)
    mesh = getbc(mesh)
    return mesh

