from timeit import default_timer as timer
from TimeIntegration import initialization, getbc, IntegrateTime
from InputVariables import nx, ngc
import matplotlib.pyplot as plt

start = timer()

def U(x,a):
    u = []
    for i in range(nx):
        u.append(x[i+1][a])
    return u

def MakeMesh(xmin,xmax,dx):
    g = [xmin]
    for i in range(nx-2):
        g.append(g[i]+dx)
    g.append(xmax)
    return g

def UpwindSolver():
    tmax, xmax, xmin, dx, x = initialization()
    t = 0
    while t < tmax:
        dt, x = IntegrateTime(x, dx)
        u = U(x, 0)
        X = MakeMesh(xmin, xmax, dx)
        plt.plot(X, u, marker='.')
        plt.show()
        t += dt
    return x, xmax, xmin, dx

UpwindSolver()

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ")
print(chrono(start,end))





