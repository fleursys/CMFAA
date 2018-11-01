import numpy as np
from InputVariables import iniCond

gamma = 1.4
nw = 3
nsnap = 30

# determines the boundary conditions, xmin, xmax, dx and tmax from the initial condition
if iniCond == 'acoustic':
    bc = 'periodic'
    xmin = 0
    xmax = 1
    tmax = 3
elif iniCond == 'shock':
    bc = 'fixed'
    xmin = -0.5
    xmax = 0.5
    tmax = 0.2

# computes the conservative variables from an array w of primitive variables (dimension nw)
def rho(w):
    return w[0]

def rhov(w):
    rv = [i * rho(w) for i in w[1:-1]]
    return rv

def e(w):
    return np.dot(w[1:-1], w[1:-1])*rho(w)/2 + w[-1]/(gamma-1)

# computes the primitive variables and the speed of sound from an array w of primitive variables (dimension nw)

def v(w):
    speed = [i/rho(w) for i in w[1:-1]]
    return speed

def p(w):
    return (gamma-1)*(-(np.dot(w[1:-1], w[1:-1])/(2*rho(w))) + w[-1])

def csound(w):
    return np.sqrt((gamma*p(w)/rho(w)))
