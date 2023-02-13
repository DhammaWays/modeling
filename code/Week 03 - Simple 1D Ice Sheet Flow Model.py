import numpy as np
import matplotlib.pyplot as plt
#from math import ??
#import math

# definition of computation parameters
n_x = 10 # grid points
domain_width = 1e6 # meters
dx = domain_width/n_x
#dz = 500
time_step = 100 # years
n_years = 2e4
n_steps = int(n_years/time_step)
flow_param = 1e4 # m horizontal / years
snow_fall = 0.5 # m/year
plot_limit = 4e3

# arrays to store calculations
elevations = np.zeros(n_x+2)
flows = np.zeros(n_x+1)

fig, ax = plt.subplots()
ax.plot(elevations)
ax.set_ylim([0,plot_limit])
plt.show(block=False)

# calculation process
for i_time in range(0, n_steps):
    for i_x in range(0, n_x+1):
        flows[i_x] = (elevations[i_x] - elevations[i_x+1])/dx*flow_param*(elevations[i_x]+elevations[i_x+1])/2/dx
        # where the first incidence of the dX variable is to
        # calculate the gradient in elevation, and the second one
        # corrects for the aspect ratio of the grid cell, how 
        # much horizontal flow translates into vertical 
        # elevation.  Be sure that you have the indexing right 
        # for your implementation of the staggered grid, and be
        # sure to get the sign right, so that flow from the first
        #cell in the list to the second would be defined as positive.  
    for i_x in range(1, n_x+1):
        # Then update the elevations by adding ( snowFall + flow[ix-1] - flow[ix] ) * timeStep.
        elevations[i_x] = elevations[i_x] + (snow_fall+flows[i_x-1]-flows[i_x])*time_step
        
    print('year', i_time*time_step)
    ax.clear()
    ax.plot(elevations)
    ax.set_ylim([0, plot_limit])
    plt.show(block=False)
    plt.pause(0.001)
    fig.canvas.draw()
ax.clear()
ax.plot(elevations)
ax.set_ylim([0, plot_limit])
plt.show()
fig.canvas.draw()
#print(flows)
#print(elevations)
