# Week 5 : Near Future Climate Model: Model for today's condition
# Autograder version

import numpy as np
#import matplotlib.pyplot as plt
import math

# Different constants

timeStep = 1 # years
eqCO2 = 280
eqCO2_ramp = 340
initCO2 = 290
CO2_exp = 0.0225
CO2RampExp = 0.01
aerosol_Wm2_now = -0.75 # -0.75 Watts/m2
watts_m2_2x = 4
climateSensitivity2x = 3 # degrees C for doubling CO2
climateSensitivityWm2 = climateSensitivity2x / watts_m2_2x
TResponseTime = 20

years = range(1900, 2100, timeStep)
bauCO2, rampCO2 = [initCO2], [initCO2]
incCO2,rfCO2,rfMask,rfCO2Ramp,rfMaskRamp,rfTot,Teq,TTrans,rfTotRamp,TeqRamp,TTransRamp = \
    [[0] for i in range(11)]

# Case: Business as usual - Humand and Industrial consumption continues!
for y in years[1:]:
    bauCO2.append(eqCO2 + (bauCO2[-1] - eqCO2) * (1 + CO2_exp * timeStep))
    incCO2.append((bauCO2[-1] - bauCO2[-2]) / timeStep)
    rfCO2.append( watts_m2_2x * math.log(bauCO2[-1]/eqCO2) / math.log(2))
    

i2015 = years.index(2018)
#aerosolCoeff = aerosol_Wm2_now / ((bauCO2[i2015] - bauCO2[i2015 -1]) / timeStep)
aerosolCoeff = aerosol_Wm2_now / incCO2[i2015]

#print(bauCO2[i2015])

for i in range(1, len(years)):
    #rfMask.append(max(incCO2[i] * aerosolCoeff, aerosol_Wm2_now))
    rfMask.append(incCO2[i] * aerosolCoeff)
    rfTot.append(rfCO2[i] + rfMask[i])
    Teq.append(rfTot[i] * climateSensitivityWm2)
    TTrans.append(TTrans[-1] + (Teq[i] - TTrans[-1]) * timeStep / TResponseTime)
 
# Rampdown

# Initial values until 2015 remain same
rampCO2 += bauCO2[1:i2015]
rfCO2Ramp += rfCO2[1:i2015]
rfMaskRamp += rfMask[1:i2015]
TTransRamp += TTrans[1:i2015]
TeqRamp += Teq[1:i2015]
rfTotRamp += rfTot[1:i2015]

for i in range(i2015, len(years)):
    rampCO2.append(rampCO2[-1] + (eqCO2_ramp - rampCO2[-1]) * (CO2RampExp * timeStep))
    rfCO2Ramp.append(watts_m2_2x * math.log(rampCO2[i] / eqCO2) / math.log(2))
    rfMaskRamp.append(0)
    rfTotRamp.append(rfCO2Ramp[i])
    TeqRamp.append(rfCO2Ramp[i] * climateSensitivityWm2)
    TTransRamp.append(TTransRamp[-1] + (TeqRamp[i] - TTransRamp[-1]) * timeStep / TResponseTime)
    
# Different plots
'''
plt.plot(years, bauCO2, '-', years, rampCO2, '-')
plt.show()
plt.plot(years, rfTot, '-', years, rfTotRamp, '-', years, rfMaskRamp, '-')
plt.show()
plt.plot(years, TTrans, '-', years, TTransRamp, '-')
plt.show()
'''

# Autograder input/output
nYear = int(input(""))
iYear = years.index(nYear)
print(bauCO2[iYear],rfCO2[iYear],Teq[iYear],TTrans[iYear])
    
    