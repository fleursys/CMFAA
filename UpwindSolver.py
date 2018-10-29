from timeit import default_timer as timer
from TimeIntegration import initialization, MakeMesh, getbc, IntegrateTime

start = timer()

def UpwindSolver():
    x = MakeMesh()
    tmax, dx, x = initialization(x)
    x = getbc(x)
    print(x)
    #t = 0
    #while t < tmax:
    #    dt, x = IntegrateTime(x,dx)
    #    t += dt
    return

UpwindSolver()

end = timer()

def chrono(start, end):
    time = end-start
    return time

print("time elapsed = ")
print(chrono(start,end))
