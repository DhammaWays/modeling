# Iterative Runaway Ice-Albedo Feedback Model

'''
Week 2 Assignment:
First we generate simple linear "parameterization" functions of planetary albedo and the latitude to
which ice forms (colder = lower latitude ice). Second, for any given value of the solar constant L,
we'll use iteration to find consistent values of albedo and T, to show the effect of the ice albedo feedback
on Earth's temperature, running away to fall into the dreaded "snowball Earth"
'''

import numpy as np
import matplotlib.pyplot as plt

#nITERS = 100
sigma = 5.67E-8 # W/m2 K4

# equation: sigma * T^4 = L/4 * (1-albedo)
def calcTemp(L, albedo, sigma=5.67E-8):    
    return pow( (L * (1 - albedo) / (4 * sigma)), 0.25 )

# Albedo and Temp linear equation: A = slope * T + Intercept
albedo_slope = -0.01
albedo_intercept = 2.8
def calcAlbedo(T, slope=albedo_slope, intercept=albedo_intercept):
    return slope * T + intercept

# Ice Lattitude and Temp linear equation: iceLat = slope * T + Intercept
iceLat_slope = 1.5
iceLat_intercept = -322.5
def calcIceLat(T, slope=iceLat_slope, intercept=iceLat_intercept):
    return slope * T + intercept

def clampMinMax(val, val_min, val_max):
    if val < val_min:
        val = val_min
    elif val > val_max:
        val = val_max
    return val

#plotType = "iter" # "L", "iter"
def runIterativeRunawayIceAlebdo(start, stop, step, nMAXITER, albedo_initial, plotType, xA, yA):
    albedo_cur = albedo_initial
    for L in range(start, stop, step):
        for it in range(nMAXITER):
            T = calcTemp(L, albedo_cur)
            albedo_cur = clampMinMax(calcAlbedo(T), 0.15, 0.65)
            iceLat = clampMinMax(calcIceLat(T), 0, 90)
            if( plotType == "iter"):
                xA.append(it)
                yA.append(T)
        # To avoid plot drawing line joining last to first
        if( plotType == "iter" ):
            xA.append(np.nan)
            yA.append(np.nan)
        
        if( plotType == "L"):
            xA.append(L)
            yA.append(T)
    # Return last calculated  
    return albedo_cur
            

# Initial Values
nITERS = 100 
albedo = 0.15 # initial guess value
albedo_l = albedo
#LRange_Min, LRange_Max, L_step = 1200, 1600, 10
LRange_Min, LRange_Max, L_step = 1150, 1350, 5

# Temperature v/s iteration
x, y = [], []
# Cooling Sweep (Sun is getting cooler): Go down
albedo = runIterativeRunawayIceAlebdo( LRange_Max, LRange_Min-1, -L_step, nITERS, albedo, "iter", x, y)
#print(albedo, x, y)
plt.plot(x, y, marker='<')
x, y = [], []
# Hot Sweep (Sun is getting hotter): Back up
albedo = runIterativeRunawayIceAlebdo( LRange_Min, LRange_Max+1, L_step, nITERS, albedo, "iter", x, y)
#print(albedo, x, y)
plt.plot(x, y, marker='>')
plt.xlabel('Iteration')
plt.ylabel('Temperature (K)')
plt.show()

# Temperature v/s L
x, y = [], []
# Cooling Sweep (Sun is getting cooler): Go down
albedo = runIterativeRunawayIceAlebdo( LRange_Max, LRange_Min-1, -L_step, nITERS, albedo, "L", x, y)
#print(albedo, x, y)
plt.plot(x, y, marker='<')
x, y = [], []
# Hot Sweep (Sun is getting hotter): Back up
albedo = runIterativeRunawayIceAlebdo( LRange_Min, LRange_Max+1, L_step, nITERS, albedo, "L", x, y)
#print(albedo, x, y)
plt.plot(x, y, marker='>')
plt.xlabel('Solar Constant (L)')
plt.ylabel('Temperature (K)')
plt.show()
