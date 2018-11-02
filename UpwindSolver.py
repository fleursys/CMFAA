from TimeIntegration import initialization, IntegrateTime
from InputVariables import var, ngc
from FixedVariables import tmax, rho, p, v, nsnap

def U(x,var, nx):
    u = []
    if var == 'density':
        for i in range(ngc, nx + ngc):
            u.append(rho(x[i]))
    elif var == 'velocity':
        for i in range(ngc, nx + ngc):
            u.append(v(x[i])[0])
    elif var == 'pressure':
        for i in range(ngc, ngc + nx):
            u.append(p(x[i]))
    else:
        print('no valid variable')
    return u

def upwind_solver(nx):
    solution = []
    time = []
    x = initialization(nx)
    t = 0
    time_snap = tmax/nsnap
    snap = time_snap
    solution.append(U(x, var, nx))
    time.append(t)
    while t < tmax:
        if t >= snap:
            dt, x = IntegrateTime(x, nx)
            t += dt
            solution.append(U(x, var, nx))
            time.append(t)
            snap += time_snap
        else:
            dt, x = IntegrateTime(x, nx)
            t += dt
    solution.append(U(x, var, nx))
    time.append(t)
    return solution, time




