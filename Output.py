from UpwindSolver import upwind_solver
from EOC import cv_and_conquer
from InputVariables import cvAnalysis, nx
from FixedVariables import xmin, xmax
from TimeIntegration import make_mesh
from timeit import default_timer as timer
import matplotlib.pyplot as plt

start = timer()

if cvAnalysis == 'false':
    solution = upwind_solver(nx)
    plot = 'solution'
elif cvAnalysis == 'true':
    convergence = cv_and_conquer()
    plot = 'convergence'

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ",chrono(start,end))

if plot == 'solution':
    X = make_mesh(xmin, xmax, nx)
    X.pop(0)
    X.pop(-1)
    for i in range(len(solution)):
        plt.plot(X, solution[i], marker='.')
        plt.show()


