import numpy 
import matplotlib.pyplot as plt
import math

nX = 10
domainWidth = 1e6 # in meters
dx = domainWidth / nX
timeStep = 100 #yrs
nYears = 20000
#nYears = float(input(''))
nSteps = int(nYears / timeStep)
flowParam = 1e4 # #horizontal / yr
snowFall = 0.5 # m/yr
plotLimit = 4000 

elevations = numpy.zeros(nX+2)
flows = numpy.zeros(nX+1)

fig, ax = plt.subplots()
ax.plot(elevations)
ax.set_ylim([0,plotLimit])
plt.show(block=False)

for iTime in range(0,nSteps):

	for ix in range(0, nX+1):
		flows[ix] = (elevations[ix] - elevations[ix+1]) / dx * flowParam * (elevations[ix] + elevations[ix+1]) /2/dx


	for ix in range(1, nX+1):
		elevations[ix] = elevations[ix] + (snowFall + flows[ix-1] - flows[ix]) * timeStep

	print( "year", iTime*timeStep)
	ax.clear()
	ax.plot(elevations)
	ax.set_ylim([0,plotLimit])
	plt.show(block = False)
	fig.canvas.draw()

print( elevations[5] )
ax.clear()
ax.plot(elevations)
ax.set_ylim([0,plotLimit])
plt.show()
fig.canvas.draw() 