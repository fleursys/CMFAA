from UpwindSolver import upwind_solver
from EOC import cv_and_conquer
from InputVariables import cvAnalysis, nx, iniCond, var
from FixedVariables import xmin, xmax
from TimeIntegration import make_mesh
from timeit import default_timer as timer
import matplotlib.pyplot as plt

start = timer()

if cvAnalysis == 'false':
    solution, time = upwind_solver(nx)
    plot = 'solution'
elif cvAnalysis == 'true':
    convergence = cv_and_conquer()
    plot = 'convergence'

end = timer()

# defines the elapsed time between the start of the computation and the end of the computation (the computational
# effort to plot the results is not counted
def chrono(start, end):
    elapsed = end-start
    return elapsed

print("time elapsed = ",chrono(start,end))

# plots the solution for every time snap
if plot == 'solution':
    X = make_mesh(xmin, xmax, nx)
    X.pop(0)
    X.pop(-1)
    for i in range(len(solution)):
        if iniCond == 'acoustic':
            plt.axis([0, 1, -0.005, 0.005])
        elif iniCond == 'fixed':
            plt.axis([-0.5, 0.5, 0, 9])
        plt.xlabel('x')
        plt.ylabel(var)
        plt.title(" t = %1.3f" %time[i])
        plt.plot(X, solution[i], marker='.')
        plt.show()


