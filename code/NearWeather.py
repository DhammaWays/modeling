import numpy as np
import matplotlib.pyplot as plt
import math

# Define the constants
dt = 5  # years
CO2_eq = 280
CO2_init = 290
CO2_exp = 0.0225
CO2_exp_ramp = 0.01
aerosol_Wm2_init = -0.75
watts_m2_2x = 4

climate_sensitivity_2x = 3
climate_sensitivity_Wm2 = climate_sensitivity_2x / watts_m2_2x
TResponseTime = 20

# initialize
years_max = 2100 # final year
years = [1900]
bauCO2 = [ CO2_init ]
incCO2 = [0]
rfCO2 = [0]
rfMask = [0]
rfCO2Ramp = [0]
rfMaskRamp = [0]
rfTot = [0]
Teq = [0]
TTrans = [0]
rampCO2 = [CO2_init]
rfTotRamp = [0]
TeqRamp = [0]
TTransRamp = [0]


# case 1: business as usual
while years[-1] < years_max:
    years.append( years[-1] + dt)
    bauCO2.append( CO2_eq + (bauCO2[-1] - CO2_eq) * (1 + CO2_exp * dt))
    incCO2.append( (bauCO2[-1]- bauCO2[-2]) / dt )
    rfCO2.append( watts_m2_2x * math.log(bauCO2[-1]/CO2_eq) / math.log(2))

i2015 = years.index(2015)
aerosol_coeff = aerosol_Wm2_init / ( (bauCO2[i2015] - bauCO2[i2015-1]) /dt )

for iY in range(1,len(years)):
    rfMask.append(max(incCO2[iY]*aerosol_coeff,aerosol_Wm2_init))
    rfTot.append(rfCO2[iY] + rfMask[iY])
    Teq.append(rfTot[iY] * climate_sensitivity_Wm2)
    TTrans.append(TTrans[-1] + (Teq[-1] - TTrans[-1])*dt/TResponseTime)

# case 2: ramp down
for iY in range(1, i2015):
    # for the years until 2015, everything stays the same
    rampCO2.append( bauCO2[iY] )
    rfCO2Ramp.append(rfCO2[iY])
    rfMaskRamp.append(rfMask[iY])
    TTransRamp.append(TTrans[iY])
    TeqRamp.append(Teq[iY])
    rfTotRamp.append(rfTot[iY])

for iY in range(i2015, len(years)):
    rampCO2.append(rampCO2[-1] + (CO2_eq*1.2 - rampCO2[-1]) * (CO2_exp_ramp * dt) )
    rfCO2Ramp.append(watts_m2_2x * math.log(rampCO2[-1]/CO2_eq) /math.log(2))
    rfMaskRamp.append(0)
    rfTotRamp.append(rfCO2Ramp[iY])
    TeqRamp.append(rfCO2Ramp[iY]*climate_sensitivity_Wm2)
    TTransRamp.append(TTransRamp[-1] + (TeqRamp[iY] - TTransRamp[-1]) * dt / TResponseTime )

# perform the plots
#plt.plot(years, bauCO2, '-', years, rampCO2, '-')
#plt.plot(years, rfTot, '-', years, rfTotRamp, '-',years, rfMaskRamp,'-')
plt.plot(years, TTrans, '-', years, TTransRamp, '-')
plt.xlabel('Time [y]')
plt.ylabel('Temperature difference [K]')
plt.show()