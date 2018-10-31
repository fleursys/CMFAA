from InputVariables import nx
import numpy as np

def AnalyticShock():
    a_sol = []
    for i in range(nx):
        a_sol.append(1)
    return a_sol

def max_error(num_sol, a_sol):
    max = 0
    for i in range(nx):
        test = np.abs(num_sol[i] - a_sol[i])
        if test > max:
            max = test
    return max






