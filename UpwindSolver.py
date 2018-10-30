from timeit import default_timer as timer
from TimeIntegration import initialization, IntegrateTime
from InputVariables import nx
from FixedVariables import xmin, xmax, dx, tmax
import matplotlib.pyplot as plt

start = timer()

def U(x,a):
    u = []
    for i in range(nx):
        u.append(x[i+1][a])
    return u

def grid(xmin,xmax,dx):
    g = [xmin]
    for i in range(nx-2):
        g.append(g[i]+dx)
    g.append(xmax)
    return g

def UpwindSolver():
    x = initialization()
    t = 0
    while t < tmax:
        dt, x = IntegrateTime(x, dx)
        t += dt
    return x

x = UpwindSolver()

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ")
print(chrono(start,end))

u = U(x, 2)
X = grid(xmin, xmax, dx)
plt.plot(X, u, marker='.')
plt.show()

