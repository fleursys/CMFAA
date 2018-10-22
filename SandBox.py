import matplotlib.pyplot as plt
import numpy as np

# define spatial dimensions
x_start = -10
x_end = + 10
delta_X = 0.1
I = int((x_end - x_start)/delta_X)

# define time discretization
t_start = 0
t_end = 10
delta_T = 0.1
N = int((t_end - t_start)/delta_T)

# define speed
v = 1

# define grid
grid = [0]*N

#define initial condition = step function
grid[0] = []
for i in range(0,int(I//2)):
    grid[0].append(-5)
for i in range(int(I//2),I):
    grid[0].append(3)

# define ghostvalues
def PeriodicBoundary_Ghosts(list, number):
    for a in range(number):
        list.insert(0, list[0])
    for a in range(number):
        list.append(list[-1])
    return list

def FTCS_advection(delta_T, delta_X, grid,v):
    GhostGrid = PeriodicBoundary_Ghosts(grid[0], 1)
    for n in range(1, N):
        grid[n] = [0] * I
        for j in range(len(GhostGrid)-2):
            grid[n][j] = GhostGrid[j+1] - v*(delta_T/(2*delta_X))*(GhostGrid[j+2]-GhostGrid[j])
        GhostGrid = PeriodicBoundary_Ghosts(grid[n], 1)
    return grid

solution = FTCS_advection(delta_T, delta_X, grid, v)

plt.plot(solution[50])
plt.show()