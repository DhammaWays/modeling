# Ice Sheet Dyanmics Model

'''
Week 3 Assignment: Autograder version

Ice flows like extra-thick molasses, downhill. The shape of the ice sheet (altitude versus distance across) is determined
by the relationship between ice surface slope and the flow rate of the ice.
'''


import numpy as np
#import matplotlib.pyplot as plt

nYears = float( input('') )

nX = 10               # number of grid points
domainWidth = 1e6      # meters
timeStep = 100         # years, 100
#nYears = 30000         # years, 50000
flowParam = 1e4        # m horizontal / yr
snowFall = 0.5         # m / y

dx = domainWidth / nX
nSteps = int(nYears / timeStep)

elevations = np.zeros(nX+2)
flows = np.zeros(nX+1)

def elevationSlope(elev, ic, dWidth=dx):
    return (elev[ic] - elev[ic+1])/dx

def cellAspectCorrection(elev, ic, dWidth=dx):
    return ((elev[ic] * 0.5 + elev[ic+1] * 0.5))/dx

# Iterate: For each step, calculate flow and update elevations
for it in range(0, nSteps):
    # Calculate flow for each cell
    for iCell in range(0, nX+1):
        flows[iCell] = elevationSlope(elevations, iCell) * flowParam * cellAspectCorrection(elevations, iCell)
    
    # Update the elevation for each cell keeping in mind eleveation is at center, zero is ghost cell
    for iCell in range(1, nX+1):
        # Elevation changes/timestep = snowfall + (inflow - outflow)
        elevations[iCell] += (snowFall + (flows[iCell-1] - flows[iCell])) * timeStep
        

# Print elevation in the center of ice sheet
print( elevations[5] )
        