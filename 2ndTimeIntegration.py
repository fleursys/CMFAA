from InputVariables import ngc
from FixedVariables import xmax, xmin, nw
from ComputeTimeStep import computeTimeStep
from TimeIntegration import getbc, get_eigenval, ConsToDiag, DiagtoCons

def get_flux(wdiag, lam, nx, dx):
    f_upwind = [[0 for l in range(nw)] for m in range(nx+ngc*2)]
    for i in range(ngc, nx + ngc):
        for j in range(nw):
            if lam[i][j] > 0:
                f_upwind[i][j] = (lam[i][j] * wdiag[i][j] - lam[i - 1][j] * wdiag[i - 1][j]) / dx
            elif lam[i][j] < 0:
                f_upwind[i][j] = (lam[i+1][j] * wdiag[i + 1][j] - lam[i][j] * wdiag[i][j]) / dx
            else:
                f_upwind[i][j] = 0
    return f_upwind

# computes the conservative variables in the next time step (after dt).
def IntegrateTime(x, nx):
    dx = (xmax - xmin) / (nx - 1)
    w = getbc(x)
    lam = get_eigenval(w, nx)
    wdiag = ConsToDiag(w, nx)
    f_upwind = get_flux(wdiag, lam, nx, dx)
    dt = computeTimeStep(lam, dx, nx)
    wdiag_adv = [[0 for l in range(nw)] for m in range(nx+ngc*2)]
    for i in range(nx + ngc*2):
        for j in range(nw):
            wdiag_adv[i][j] = wdiag[i][j] - dt*f_upwind[i][j]
    w_adv = DiagtoCons(w, wdiag_adv, nx)
    return dt, w_adv