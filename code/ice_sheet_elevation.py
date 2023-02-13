# Ice Sheet Dyanmics Model

'''
Week 3 Assignment

Ice flows like extra-thick molasses, downhill. The shape of the ice sheet (altitude versus distance across) is determined
by the relationship between ice surface slope and the flow rate of the ice.

Flow = dE/dX * Constant

Ice flows, like any other fluid, only very slowly, due to its high viscosity. The force driving the flow is differences
in pressure in the interior of the ice, which arise from differences in the elevation of the ice surface. If you make
a pile of ice, it's like a pile of molasses, in that it will flow outward and flatten itself out. The model is formulated
in one dimension, on a horizontal grid. Start with 10 grid cells. Let them span a horizontal distance of 1000 km,
or 10^6 meters. Each grid cell will have an elevation of ice. Flow between adjacent cells depends on the difference
between their elevations. Snow falls equally on all grid cells. Also, vitally important to this type of problem, 
is a boundary condition. We'll assume that the ice sheet is confined to a landmass like Greenland or Antarctica,
so that the thickness of the ice at the boundaries has to be zero. Ice flows into the ocean and disappears, both 
in reality and in our model. 

The model will step forward in time, using a time step of 100 years. It will begin to rise uniformly, but the elevations
at the edges will be eroded by flow to the ocean. Eventually it will reach a steady state, where snowfall in each grid cell
is balanced by flow from the grid cell, which is all determined by the slope of the ice surface. 

'''

import numpy as np
import matplotlib.pyplot as plt

nX = 20                # number of grid points
domainWidth = 1e6      # meters
timeStep = 40         # years, 100
nYears = 20000         # years, 50000
flowParam = 1e4        # m horizontal / yr
snowFall = 0.5         # m / y

dx = domainWidth / nX
nSteps = nYears // timeStep

elevations = np.zeros(nX+2)
flows = np.zeros(nX+1)

def elevationSlope(elev, ic, dWidth=dx):
    return (elev[ic] - elev[ic+1])/dx

def cellAspectCorrection(elev, ic, dWidth=dx):
    return ((elev[ic] * 0.5 + elev[ic+1] * 0.5))/dx

# Setup realtime plot
fig, ax = plt.subplots()
ax.plot(elevations)
ax.set_ylim([0,4000])
plt.show(block=False)

# Iterate: For each step, calculate flow and update elevations
for it in range(0, nSteps):
    # Calculate flow for each cell
    for iCell in range(0, nX+1):
        flows[iCell] = elevationSlope(elevations, iCell) * flowParam * cellAspectCorrection(elevations, iCell)
    
    # Update the elevation for each cell keeping in mind eleveation is at center, zero is ghost cell
    for iCell in range(1, nX+1):
        # Elevation changes/timestep = snowfall + (inflow - outflow)
        elevations[iCell] += (snowFall + (flows[iCell-1] - flows[iCell])) * timeStep
     
    # Keep updating plot
    print("Year: ", it * timeStep)
    ax.clear()
    ax.plot( elevations )
    plt.show( block=False )
    plt.pause(0.001)
    fig.canvas.draw()

# Draw final elevation plot

plt.plot(elevations, marker='^')
plt.ylabel('Snow Elevation (m)')
plt.xlabel('Cell (width={} m)'.format(dx))
plt.show()
        