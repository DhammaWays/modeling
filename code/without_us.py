# Week 5 : Near Future Climate Model: Model for today's condition

'''
Description of the Model Formulation
Simulate the near-term future of Earth’s climate as a function of the three main uncertainties: 
-future CO2 emission, 
-the climate sensitivity, 
-and the current cooling impact of industrial aerosols and short-lived greenhouse gases (the masking effect). 

The simulation has two stages, one with people releasing CO2 and generating aerosols, which we'll call Business-as-Usual,
and a second, diverging from the first, in which we suddenly stop both of those things, at which point the CO2 largely
persists, while the aerosols go away immediately. 

Numerical Guts of the Model
1. In either framework, you will need a list of years, from 1900 to 2100 in steps of 1 year. 
2. Business-as-usual CO2 is increasing exponentially, following a function of the form
pCO2 (time 2) = 280 + (pCO2 (time 1)-280) * (1 + A * D_time)
where D_time is the time step, time_1 and time_2 are time points (start from a known pCO2 at time_1 and calculate what it
will be at time_2). A is a tunable growth rate parameter equal to a value of 0.0225 / year, gives an atmospheric rise rate
of 2.5 ppm when the pCO2 value is 400 ppm; a reasonable fit to the current pCO2 value and rate of rise. 
3. For each pCO2 level, calculate the radiative forcing using the formula
RF from CO2 = 4 [Watts/m2] * ln( pCO2 / 280 ) / ln ( 2 ), where the last factor, ln(pCO2/280) / ln(2), gives the number
of doublings of the CO2 concentration over an initial value of280. 
4. We need to account for the radiative impact of short-lived industrial stuff, primarily the sum of cooling from sulfate aerosols and warming from short-lived greenhouse gases like methane. But the difficulty is that we don't know very well what that number is. Assume that the intensity of the masking (how strong it is in Watts / m2) depends on the rate of industrial activity, which we can estimate from the rate of pCO2 rise. In other words, use a function of the form
RF masking scaled = B * (pCO2 (year 2015) – pCO2 (year 2014))/( 1 year), where you can get the value of B by fitting to
a radiative forcing estimate for the present day. In other words, to get a masking effect of –0.75 Watts / m2, when the
pCO2 rise is 2.5 ppm / year, B would have to be –0.3. Solve for B in the code based on the assumption that the present-day
masking RF is -0.75 Watts/m2. Later you'll change that assumption, so the code should know how to calculate B on its own. 
The model looks more reasonable in the future if we assume that this function only works for the past, but that further 
increases in CO2 emissions don't lead to still further increases in aerosols, as in 
RF masking = max( RF masking scaled, RF masking in 2015 )
It is easy for coal plants to clean up their sulfur, which leads to the aerosols (and lots of air quality problems), without cutting the carbon emissions. 
5. Compute the total radiative forcing as the sum from CO2 and from masking. 
6. Compute the equilibrium temperature change that you would get from that RF, using a climate sensitivity Delta_T_2x of
3 degrees per doubling of CO2, but to get that number into units of Watts/m2 rather than CO2 doublings (so it can deal
with the effect of aerosols), use the fact that doubling CO2 gives 4 Watts/m2 of forcing, so that
climate_sensitivity_Watts_m2 = climate_sensitivity_2x / 4 [Watts/m2 per doubling CO2]
7. Compute the evolution of temperature by relaxing toward equilibrium on a time scale of 100 years, in other words,taking
change in T per timestep = (T(equilibrium) – T) / t_response_time (20 years) * X years / timestep

That was business as usual. Now for The World Without Us. 
8. Beginning in the year when pCO2 hits 400, make a new column or list with pCO2 values that decrease rather than increase,
following a somewhat fudged form
change in CO2 per timestep = ( 340 - CO2 ) * ( 0.01 * timestep ), where the CO2 concentration is now relaxing toward a
higher concentration than the initial, due to the long-term change in ocean chemistry (acidification). The time scale for
the CO2 invasion into the ocean is 100 years (1 / 0.01), which is a composite of fast equilibration times with the surface
ocean and the land biosphere, and slower with the deep ocean. 
9. Compute the RF from CO2 from this as you did for the last pCO2 column. 
10. Compute the equilibrium temperature from the RF from CO2. There is no masking RF for this stage, because that all goes
away quickly. 
11. Compute the time-evolving temperature from this period from the equilibrium temperature as you did before. 

'''

import numpy as np
import matplotlib.pyplot as plt
import math

# Different constants
timeStep = 1 # years
eqCO2 = 280
eqCO2_ramp = 340
initCO2 = 290
CO2_exp = 0.0225
CO2RampExp = 0.01
aerosol_Wm2_now = -0.75 # W/m2
watts_m2_2x = 4
climateSensitivity2x = 3 # degrees C for doubling CO2
climateSensitivityWm2 = climateSensitivity2x / watts_m2_2x
TResponseTime = 20

# Initial conditions
years = range(1900, 2100, timeStep)
bauCO2, rampCO2 = [ initCO2 ], [initCO2]
incCO2,rfCO2,rfMask,rfCO2Ramp,rfMaskRamp,rfTot,Teq,TTrans,rfTotRamp,TeqRamp,TTransRamp = \
    [[0] for i in range(11)]

# Case: Business as usual - Humand and Industrial consumption continues!
for y in years[1:]:
    bauCO2.append(eqCO2 + (bauCO2[-1] - eqCO2) * (1 + CO2_exp * timeStep))
    incCO2.append((bauCO2[-1] - bauCO2[-2]) / timeStep)
    rfCO2.append( watts_m2_2x * math.log(bauCO2[-1]/eqCO2) / math.log(2))
    
i2015 = years.index(2015)
aerosolCoeff = aerosol_Wm2_now / ((bauCO2[i2015] - bauCO2[i2015 -1]) / timeStep)

for i in range(1, len(years)):
    rfMask.append(max(incCO2[i] * aerosolCoeff, aerosol_Wm2_now))
    rfTot.append(rfCO2[i] + rfMask[i])
    Teq.append(rfTot[i] * climateSensitivityWm2)
    TTrans.append(TTrans[-1] + (Teq[i] - TTrans[-1]) * timeStep / TResponseTime)

# Case: Rampdown: World without us!

# Initial values just before 2015 remain same
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
    
# To answer quiz questions
#print(climateSensitivity2x, TTrans[i2015])
#TCompare = [(y+2014, t1, t2) for y, (t1, t2) in enumerate(zip(TTrans[i2015-1:i2015+40], TTransRamp[i2015-1:i2015+40]))]
#print(TCompare)

# Different plots
plt.plot(years, bauCO2, '-', years, rampCO2, '-')
plt.show()
plt.plot(years, rfTot, '-', years, rfTotRamp, '-', years, rfMaskRamp, '-')
plt.show()
plt.plot(years, TTrans, '-', years, TTransRamp, '-')
plt.show()

        