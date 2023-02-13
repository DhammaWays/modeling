# Week 4: 2D Gridded Shallow Water Model: Autograder version

"""

The first section of the code contains setup and initialization 
information.  Leave it alone for now, and you can play with them later 
after you get the code filled in and running without bugs.  

"""

# Set up python environment.  numpy and matplotlib will have to be installed
# with the python installation.

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.ticker as tkr
import math

# Input for Autograder
ncol = 3
nSlices, iRowOut, iColOut = input("").split()
nSlices, iRowOut, iColOut = [ int(nSlices), int(iRowOut), int(iColOut) ]         
ntAnim = 1          
horizontalWrap = True 
interpolateRotation = False  # or True, either works
textOutput = False
plotOutput = False

dT = 600    # seconds
G = 9.8e-4  # artificially low to allow a long time step

# Grid and Variable Initialization -- stuff you might play around with

#ncol = 10         # grid size (number of cells)
nrow = ncol
#nSlices = 20         # 2, maximum number of frames to show in the plot
#ntAnim = 10          # 1, number of time steps for each frame

#horizontalWrap = False # True, determines whether the flow wraps around, connecting
                       # the left and right-hand sides of the grid, or whether
                       # there's a wall there. 
#interpolateRotation = True
rotationScheme = "PlusMinus"   # "WithLatitude", "PlusMinus", "Uniform"

# Note: the rotation rate gradient is more intense than the real world, so that
# the model can equilibrate quickly.

windScheme = ""  # "Curled", "Uniform"
initialPerturbation = "Tower"    # "Tower", "NSGradient", "EWGradient"
#textOutput = False # False
#plotOutput = True
arrowScale = 30

#dT = 600 # seconds
#G = 9.8e-4 # m/s2, hacked (low-G) to make it run faster
HBackground = 4000 # meters

dX = 10.E3 # 10.E3 meters, small enough to respond quickly.  This is a very small ocean
# on a very small, low-G planet.  

dxDegrees = dX / 110.e3
flowConst = G  # 1/s2
dragConst = 1.E-6  # about 10 days decay time
meanLatitude = 30 # degrees

# Here's stuff you probably won't need to change

latitude = []
rotConst = []
windU = []
for irow in range(0,nrow):
    if rotationScheme == "WithLatitude":
        latitude.append( meanLatitude + (irow - nrow/2) * dxDegrees )
        rotConst.append( -7.e-5 * math.sin(math.radians(latitude[-1]))) # s-1
    elif rotationScheme == "PlusMinus":
        rotConst.append( -3.5e-5 * (1. - 0.8 * ( irow - (nrow-1)/2 ) / nrow )) # rot 50% +-
    elif rotationScheme == "Uniform":
        rotConst.append( -3.5e-5 ) 
    else:
        rotConst.append( 0 )

    if windScheme == "Curled":
        windU.append( 1e-8 * math.sin( (irow+0.5)/nrow * 2 * 3.14 ) ) 
    elif windScheme == "Uniform":
        windU.append( 1.e-8 )
    else:
        windU.append( 0 )
itGlobal = 0

# Add a constant for ghost row
rotConstV = np.concatenate([[1], rotConst]) 
windU = np.concatenate([[0], windU])

# Initialize velocities, elevation, rotation and needed derivatives
nrow_g = nrow + 1 # Ghost row added
ncol_g = ncol + 1 # Ghost column added

U = np.zeros((nrow_g, ncol_g))
V = np.zeros((nrow_g, ncol_g))
H = np.zeros((nrow+1, ncol+1))
dUdT = np.zeros((nrow_g, ncol_g))
dVdT = np.zeros((nrow_g, ncol_g))
dHdT = np.zeros((nrow_g, ncol_g))
dHdX = np.zeros((nrow_g, ncol_g))
dHdY = np.zeros((nrow_g, ncol_g))
dUdX = np.zeros((nrow_g, ncol_g))
dVdY = np.zeros((nrow_g, ncol_g))
rotV = np.zeros((nrow_g,ncol_g)) # interpolated to u locations
rotU = np.zeros((nrow_g,ncol_g)) #              to v

#print(V)
#print(rotV)
#print(V[1:,:])
    
midCell = int(ncol/2)
if initialPerturbation == "Tower":
    H[midCell+1,midCell] = 1
elif initialPerturbation == "NSGradient":
    H[0:midCell+1,:] = 0.1
elif initialPerturbation == "EWGradient":
    H[:,0:midCell] = 0.1

"""
This is the work-horse subroutine.  It steps forward in time, taking ntAnim steps of
duration dT.  
"""
    
def animStep():    

    global stepDump, itGlobal
    global U, V, H, dUdT, dVdT, dHdT, dHdX, dUdX, dVdY, rotV, rotU

    # Time Loop
    for it in range(0,ntAnim):
    # Here is where you need to build some code
    
        # Encode Longitudinal Derivatives Here
        dHdX = (H - np.roll(H, 1, axis=1))/dX  # H[i] - H[i-1]
        dUdX = (np.roll(U, -1, axis=1) - U)/dX # U[i+1] - U[i]
    
        # Encode Latitudinal Derivatives Here
        dHdY = (H - np.roll(H, 1, axis=0))/dX  # H[j] - H[j-1]
        dVdY = (np.roll(V, -1, axis=0) - V)/dX # V[j+1] - V[j]
         
        # Calculate the Rotational Terms Here
        rotU = U * rotConstV[:, None]  # Used for V
        rotV = V * rotConstV[:, None] # Used for U

        # Assemble the Time Derivatives Here
        dUdT = rotV - flowConst * dHdX - dragConst * U + windU[:, None]
        dVdT = -rotU - flowConst * dHdY - dragConst * V
        dHdT = -(dUdX + dVdY) * (HBackground/dX)

        '''
        print("step:", itGlobal)
        print("dHdX:")
        print(dHdX)
        print("dUdX:")
        print(dUdX)

        print("dHdY:")
        print(dHdY)
        print("dVdY:")
        print(dVdY)
        
        print("rotU")
        print(rotU)
        print("rotV")
        print(rotV)
        
        print("dUDT:")
        print(dUdT)
        
        print("Before U:")
        print(U)
 
        print("dVDT:")
        print(dVdT)
 
        print("Before V:")
        print(V)
        '''
       
        # Step Forward One Time Step
        U += dUdT * dT
        V += dVdT * dT
        H += dHdT * dT
        
        #print("After U:")
        #print(U)
        #print("After V:")
        #print(V)
    
        # Update the Boundary and Ghost Cells
        U[0], V[0] = 0., 0. # Velocities at north wall are zero
        if( horizontalWrap ):
            U[:,ncol] = U[:,0]
            H[:,ncol] = H[:,0]
        else:
            U[:,0], U[:,ncol] = 0., 0. # Velocities at east and west are zero

#   Now you're done

    itGlobal = itGlobal + ntAnim

def firstFrame():
    global fig, ax, hPlot
    fig, ax = plt.subplots()
    ax.set_title("H")   
    #hh = H[:,0:ncol]
    hh = H[1:,0:ncol]
    loc = tkr.IndexLocator(base=1, offset=1)
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    grid = ax.grid(which='major', axis='both', linestyle='-')
    hPlot = ax.imshow(hh, interpolation='nearest', clim=(-0.5,0.5))   
    plotArrows()
    #plt.show(block=False) 
    plt.show(block=True) 

def plotArrows():
    global quiv, quiv2
    xx = []
    yy = []
    uu = []
    vv = []
    Ur = U[1:]
    Vr = V[:,:-1]
    for irow in range( 0, nrow ):
        for icol in range( 0, ncol ):
            xx.append(icol - 0.5)
            yy.append(irow )
            #uu.append( U[irow,icol] * arrowScale )
            uu.append( Ur[irow,icol] * arrowScale )
            vv.append( 0 )
    quiv = ax.quiver( xx, yy, uu, vv, color='white', scale=1)
    for irow in range( 0, nrow ):
        for icol in range( 0, ncol ):
            xx.append(icol)
            yy.append(irow - 0.5)
            uu.append( 0 )
            #vv.append( -V[irow,icol] * arrowScale )
            vv.append( -Vr[irow,icol] * arrowScale )

    quiv2 = ax.quiver( xx, yy, uu, vv, color='white', scale=1)

def updateFrame():
    global fig, ax, hPlot, quiv, quiv2
    #hh = H[:,0:ncol]
    hh = H[1:,0:ncol]
    hPlot.set_array(hh)
    #quiv.remove()    
    #quiv2.remove()
    #plotArrows()
    #plt.pause(0.00001)
    #fig.canvas.draw()
    plt.show()
    print("Time: ", math.floor( itGlobal * dT / 86400.*10)/10, "days")

def textDump():
    print("time step ", itGlobal)    
    print("H", H)
    print("dHdX" )
    print( dHdX)
    print("dHdY" )
    print( dHdY)
    print("U" )
    print( U)
    print("dUdX" )
    print( dUdX)
    print("rotV" )
    print( rotV)
    print("V" )
    print( V)
    print("dVdY" )
    print( dVdY)
    print("rotU" )
    print( rotU)
    print("dHdT" )
    print( dHdT)
    print("dUdT" )
    print( dUdT)
    print("dVdT" )
    print( dVdT)

if textOutput is True:
    textDump()
if plotOutput is True:
    firstFrame()
for i_anim_step in range(0,nSlices):
    animStep()
    if textOutput is True:
        textDump()
    if plotOutput is True:
        #updateFrame()
        firstFrame()
#print(dHdT[0,1])

# Output for autograder
#print("H:", H)
#print("dHdT:", dHdT)
#print("U:", U)
#print("V:", V)
#print("rotU:", rotU)
iRowOut += 1 # Top one is ghost row!
print(H[iRowOut,iColOut],dHdT[iRowOut,iColOut],U[iRowOut,iColOut],V[iRowOut,iColOut],rotU[iRowOut,iColOut])
