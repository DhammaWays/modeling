# Naked Planet Model
'''
Week 1 Assignment: Naked Planet Model for temperature change

The temperature of a planet is determined by balancing energy fluxes into and out of a planet. Incoming solar heat
is determined by:
L * (1-albedo) / 4, and outgoing infrared is calculated as epsilon * sigma * T^4.

The goal is to numerically simulate how the planetary temperature of a naked planet would change through time as it
approaches equilibrium (the state at which it stops changing, which we calculated before). The planet starts with some
initial temperature. The “heat capacity” (units of Joules / m2 K) of the planet is set by a layer of water which absorbs
heat and changes its temperature. If the layer is very thick, it takes a lot more heat (Joules) to change the temperature.
The differential equation you are going to solve is: 
dHeatContent/dt = L*(1-alpha)/4 - epsilon * sigma * T^4, where the heat content is related to the temperature by the heat capacity 

T[K] = HeatContent [J/m2] / HeatCapacity [J/m2 K]
The numerical method is to take time steps, extrapolating the heat content from one step to the next using the incoming
and outgoing heat fluxes, same as you would balance a bank account by adding all the income and subtracting all the
expenditures over some time interval like a month. The heat content "jumps" from the value at the beginning of the
time step, to the value at the end, by following the equation:
HeatContent(t+1) = HeatContent(t) + dHeatContent/dT * TimeStep

This scheme only works if the time step is short enough that nothing too huge happens through the course of the time step. 
Set the model up with the following constants:
7
timeStep = 100           # years
waterDepth = 4000        # meters
L = 1350                 # Watts/m2
albedo = 0.3
epsilon = 1
sigma = 5.67E-8          # W/m2 K4

'''

import numpy as np
#import matplotlib.pyplot as plt

nSTEPS = int(input(""))
#nSTEPS = 100
timeStep = 100 # years
waterDepth = 4000 # meters
L = 1350 # Watts/m2
albedo = 0.3
epsilon = 1
sigma = 5.67E-8 # W/m2 K4
secondsPerYear = 3.14e7 # 365*24*60*60

heatCapacity = waterDepth * 4.2E6 # J/K m2
TK = [0.]
heatContent = heatCapacity * TK[0]
heatIn = L * (1- albedo)/4 # remains constant
heatOut = 0
timeYears = np.arange(0, (nSTEPS+1)*timeStep, timeStep)
#print(heatCapacity)
for it in range(0, nSTEPS):
    heatOut = epsilon * sigma * pow(TK[-1], 4)
    #print(timeYears[it], heatOut)
    heatContent += (heatIn - heatOut) * timeStep * secondsPerYear
    TK.append(heatContent / heatCapacity)
    #print(timeYears[it+1], TK[-1], heatContent, heatIn, epsilon * sigma * pow(TK[-1], 4))

# Last temperature (K) and heat flux out
heatOut = epsilon * sigma * pow(TK[-1], 4)
print(TK[-1], heatOut)    
# plot(x=timeYears, y=TK)
#plt.plot(timeYears, TK)
#plt.xlabel("Time in Years")
#plt.ylabel("Temperature")
#plt.show()()