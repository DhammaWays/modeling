# Week 4: 2D Gridded Shallow Water Model: 

"""
Geostrophic Drifting Rossby Wave Simulation

The flow in the model is driven by differences in the elevation of the water surface, which arise from the flow itself.
The flow is also altered by rotation (as on a rotating planet), potentially by wind, and by friction. The water is assumed
to be homogeneous in the vertical, no differences in temperature, density, or velocity are tracked. Because of the vertical
homogeneity, these are called the "shallow water equations".

In its simplest configuration, this model can start with an initial "hill" of high water level in the center of the
computational grid. Water starts to flow outward from the hill, but it rotates to the side, and so tends to find a
pattern where the flow is going around and around the hill, rather than simply flowing straight down as it would if
there were no rotation.

Note that steps have been taken in the template file to accelerate the evolution of the model, so that the relatively
slow Python language can make plots evolve at a visually satisfying fast rate. These adjustments are (1) low gravity,
which slows down the waves and allows a longer time step, (2) a small ocean with grid cells that are not too large,
so that the amount of water in them can change quickly, and (3) a much steeper change in rotation rate with latitude,
as if the ocean was on a very small planet. 

For the simple 3x3 case, the placement of the velocities and elevations is as shown below. Longitudinal velocities (u)
are defined on the left-hand cell faces, and latitudinal velocities (V) on the cell tops. The elevation of the fluid in
each box (H) is defined in the cell centers.

----V(00)-------V(01)--------V(02)----
  |           |            |           |
U(00) H(00) U(01) H(01)  U(02) H(02) [U(03)]
  |           |            |           |
  ----V(10)-------V(11)--------V(12)----
  |           |            |           |
U(10) H(10) U(11) H(11)  U(12) H(12) [U(13)]
  |           |            |           |  
  ----V(20)-------V(21)--------V(22)----
  |           |            |           |
U(20) H(20) U(21) H(21)  U(22) H(22) [U(23)]
  |           |            |           |  
  ---[V(30)]-----[V(31)]------[V(32)]---
  
Variables that are enclosed with square brackets in this diagram are "ghost" variables. They aren't computed as part
of the real grid, but are used to make it simpler to calculate differences, say between the North/South velocities (V)
at the top and the bottom of each cell.

The Differential Equations
dU/dT = C_rotation * V - C_flow * dH/dX - C_drag * U + C_wind
dV/dT = - C_rotation * U - C_flow * dH/dY - C_drag * V
dH/dT = - ( dU/dX + dV/dY ) * Hbackground / dX

The terms like dU/dT denote derivatives, in this case the rate of change of the velocity (U) with time (t). C_ terms are
constants, for rotation, induction of flow, drag, or wind input. Notice how the rotation term transfers energy between the
U and V velocities. Terms like dV/dY or dH/dX are spatial derivatives, how much the velocity (V) changes with latitude (y),
for example. Ultimately these "time tendency" (e.g. dU/dT) terms will be used to update the values according to, for example

U(time+1) = U(time) + dU/dT * delta_t, where delta_t is a time step. 

The first equation can be interpreted as a list of things (on the right hand side) which tend to make the velocity (u)
change with time (dU/dT). Rotation changes the flow direction, moving velocity between the two directions U and V according
to a rotational constant, C_rotate, where a higher number (away from zero) would rotate faster, and the sign of the constant determines the direction of rotation. The term dh/dX tells whether the sea surface is sloping; if it is, it drives the flow to accellerate according to a flow constant C_flow.
Drag slows the flow down, the larger the flow (U), the faster the slowdown (dU/dT). And finally we'll set up to allow wind
to blow, in an East/West (U) direction, driving circulation.

Mapping the Equations Numerically onto the Grid
The differential equations, above, apply to a continuous fluid, where you can imagine a slope of the sea surface as a
tangent to a wavy surface. Numerically, we cast these equations onto much coarser systems of boxes, and derive our pressure
driver, for example, from differences of heights between adjacent boxes. We need to do this paying attention to how the
variables (U, V, and H) are arrayed in space on the grid (above).

For example, in the equation for dU/dT, above, the flow is driven by a sloping sea surface (dH/dX). Looking at the grid 
diagram, the slope in the sea surface (H), appropriate to the second U in from the left, U(01), would span the position 
of U(01), to be H(01) - H(00), divided by the grid spacing delta_x. 

  |           |            |           |
U(00) H(00) U(01) H(01)  U(02) H(02) [U(03)]
  |           |            |           |

In the code, we want to have grids of the three important variables U, V, and H, and also we'll construct arrays of intermediate variables like dH/dX and dV/dY. They will be indexed in the grid, as in
 dU/dT[irow,icol] = flowConst * dH/dX[irow,icol] + ...

Rotation
The effect of Earth's rotation is to transfer velocity between the U and V directions. There are two ways to calculate
this effect: an easy but somewhat biased way which will mostly work but lead to some weird flow patterns, and a better way.
You don't have to do either one to pass the Code Check, but there are optional code checkers for both schemes if you want to
try your hand. The two formulations diverge in longer simulations, where the simpler method will generate diagonal "stripes"
of water level (waves) across the grid, as an artifact of its less-accurate method.

Easier First Method: Loop over all rows and columns, and calculate: 
rotU[irow,icol] = rotConst[irow] * U[irow,icol]

This array will be used to update V[irow,icol]. rotConst is a function of latitude so that rotation can be stronger in
high latitudes, like rotation on a real sphere.


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import math

# Grid and Variable Initialization -- stuff you might play around with

ncol = 10         # grid size (number of cells)
nrow = ncol
nSlices = 20         # 2, maximum number of frames to show in the plot
ntAnim = 2          # 1, number of time steps for each frame

horizontalWrap = False # True, determines whether the flow wraps around, connecting
                       # the left and right-hand sides of the grid, or whether
                       # there's a wall there. 
interpolateRotation = True
rotationScheme = "PlusMinus"   # "WithLatitude", "PlusMinus", "Uniform"

# Note: the rotation rate gradient is more intense than the real world, so that
# the model can equilibrate quickly.

windScheme = ""  # "Curled", "Uniform"
initialPerturbation = "Tower"    # "Tower", "NSGradient", "EWGradient"
textOutput = False # False
plotOutput = True
arrowScale = 30

dT = 600 # seconds
G = 9.8e-4 # m/s2, hacked (low-G) to make it run faster
HBackground = 4000 # meters

dX = 1.E3 # 10.E3 meters, small enough to respond quickly.  This is a very small ocean
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
# Adding a top ghost row and last ghost column to make calculations easier!
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

# Pertub initial hieight a bit to get the simulation going!    
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
    
        # Taking advantage of numpy vectorize operations for faster code (and easier to read) rather slow for loop!
    
        # Encode Longitudinal Derivatives Here
        dHdX = (H - np.roll(H, 1, axis=1))/dX  # H[i] - H[i-1]
        dUdX = (np.roll(U, -1, axis=1) - U)/dX # U[i+1] - U[i]
    
        # Encode Latitudinal Derivatives Here
        dHdY = (H - np.roll(H, 1, axis=0))/dX  # H[j] - H[j-1]
        dVdY = (np.roll(V, -1, axis=0) - V)/dX # V[j+1] - V[j]
         
        # Calculate the Rotational Terms Here
        rotU = U * rotConstV[:, None]  # Used for V
        rotV = V * rotConstV[:, None]  # Used for U

        # Assemble the Time Derivatives Here
        dUdT = rotV - flowConst * dHdX - dragConst * U + windU[:, None]
        dVdT = -rotU - flowConst * dHdY - dragConst * V
        dHdT = -(dUdX + dVdY) * (HBackground/dX)
       
        # Step Forward One Time Step
        U += dUdT * dT
        V += dVdT * dT
        H += dHdT * dT
    
        # Update the Boundary and Ghost Cells
        U[0], V[0] = 0., 0. # Velocities at north wall are zero
        if( horizontalWrap ): # Flow wraps around from leftmost to rightmost
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
    hh = H[1:,0:ncol] # Top row is ghost row!
    loc = tkr.IndexLocator(base=1, offset=1)
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    grid = ax.grid(which='major', axis='both', linestyle='-')
    hPlot = ax.imshow(hh, interpolation='nearest', clim=(-0.5,0.5))   
    plotArrows()
    plt.show(block=False)  

def plotArrows():
    global quiv, quiv2
    xx = []
    yy = []
    uu = []
    vv = []
    Ur = U[1:]     # Top row is ghost row
    Vr = V[:,:-1]  # Righmost column is ghost column
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
    hh = H[1:,0:ncol] # Top row is a ghost row!
    hPlot.set_array(hh)
    quiv.remove()    
    quiv2.remove()
    plotArrows()
    plt.pause(0.00001)
    fig.canvas.draw()
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
textDump()

