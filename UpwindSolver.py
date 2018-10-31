from timeit import default_timer as timer
from TimeIntegration import initialization, IntegrateTime, make_mesh
from InputVariables import nx, var, ngc, cvAnalysis
from FixedVariables import xmin, xmax, dx, tmax, rho, p, v, nsnap
from EOC import AnalyticShock, max_error
import matplotlib.pyplot as plt

start = timer()

def U(x,var):
    u = []
    if var == 'density':
        for i in range(ngc, nx + ngc):
            u.append(rho(x[i]))
    elif var == 'velocity':
        for i in range(nx):
            u.append(v(x[i])[0])
    elif var == 'pressure':
        for i in range(ngc, ngc + nx):
            u.append(p(x[i]))
    else:
        print('no valid variable')
    return u

def UpwindSolver():
    solution = []
    x = initialization()
    t = 0
    time_snap = tmax/(nsnap)
    solution.append(U(x, var))
    while t < tmax:
        if t >= time_snap:
            dt, x = IntegrateTime(x, dx)
            t += dt
            solution.append(U(x, var))
            time_snap += time_snap
        else:
            dt, x = IntegrateTime(x, dx)
            t += dt
    solution.append(U(x, var))
    return solution

def cv_and_conquer():
    num_sol = UpwindSolver()[-1]
    a_sol = AnalyticShock()
    err_1 = max_error(num_sol, a_sol)
    err_2 = 3
    return err_1 - err_2


Solution = UpwindSolver()

if cvAnalysis == 'true':
    a = 2

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ")
print(chrono(start,end))

