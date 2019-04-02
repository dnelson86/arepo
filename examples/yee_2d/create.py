""" @package ./examples/Yee_2d/create.py
Code that creates 2d Yee vortex initial conditions
parameters are identical to Pakmor et al. (2016), MNRAS 455, 1134

created by Rainer Weinberger, last modified 21.02.2019 -- comments welcome
"""

#### load libraries
import sys    # system specific calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/Yee_2d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(10.0)
if len(sys.argv) > 3:
  CellsPerDimension = IntType(sys.argv[3])
else:
  CellsPerDimension = IntType(50)

#### parameters
Tinf = FloatType(1.0)
beta = FloatType(5.0)
gamma = FloatType(1.4)  ## note: this has to be consistent with the parameter settings for Arepo!

""" set up grid: equidistant polar 2d grid; similar to Pakmor et al (2016) """
d_ring = Boxsize / CellsPerDimension
## place central cell
Radius = [0.0]
xPosFromCenter = [0.0]
yPosFromCenter = [0.0]

i_ring = 0
while d_ring * FloatType(i_ring) <= 0.71 * Boxsize:
    i_ring += 1
    ## place i_ring cells at this distance
    phi = np.random.uniform(0,2.0*np.pi) ## random starting angle
    #print "i_ring", i_ring
    for i_cell in np.arange( IntType(2.0*np.pi*i_ring) ):
        radius_this_cell = d_ring * FloatType(i_ring)
        #print radius_this_cell
        xcoord = radius_this_cell * np.sin(phi)
        ycoord = radius_this_cell * np.cos(phi)
        
        ## only include the cell if coordinates are within the box
        if xcoord >= -0.5*Boxsize and xcoord < 0.5*Boxsize and ycoord >= -0.5*Boxsize and ycoord < 0.5*Boxsize:
            Radius.append( radius_this_cell )
            xPosFromCenter.append( xcoord )
            yPosFromCenter.append( ycoord )
            
        phi += (2.0*np.pi*i_ring) / IntType(2.0*np.pi*i_ring) / FloatType(i_ring)

## convert to numpy arrays now that number of cells is known
Radius = np.array(Radius, dtype=FloatType)
xPosFromCenter = np.array(xPosFromCenter, dtype=FloatType)
yPosFromCenter = np.array(yPosFromCenter, dtype=FloatType)
NumberOfCells = Radius.shape[0]

## set up structure for positions (in code coordinates, i.e. from 0 to Boxsize)
Pos = np.zeros([NumberOfCells, 3])
Pos[:,0] = xPosFromCenter + 0.5 * Boxsize
Pos[:,1] = yPosFromCenter + 0.5 * Boxsize

""" set up hydrodynamical quantitites """
Temperature = Tinf - ( ( gamma - 1.0 ) * beta * beta / ( 8.0 * np.pi * np.pi * gamma ) * np.exp( 1.0 - Radius*Radius ) )
exponent = 1.0 / ( gamma - 1.0 )

Density = Temperature**exponent
Vel = np.zeros([NumberOfCells, 3], dtype=FloatType)
Vel[:,0] = -0.5 * yPosFromCenter * beta / np.pi * np.exp( 0.5 * (1.0 - Radius * Radius) )
Vel[:,1] = 0.5 * xPosFromCenter * beta / np.pi * np.exp( 0.5 * (1.0 - Radius * Radius) )
Uthermal = Temperature / (gamma - 1.0)

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(simulation_directory+'/IC.hdf5', 'w')

## create hdf5 groups
header = IC.create_group("Header")
part0 = IC.create_group("PartType0")

## header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype = IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype = IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype = IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## copy datasets
part0.create_dataset("ParticleIDs", data = np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data = Pos)
part0.create_dataset("Masses", data = Density)
part0.create_dataset("Velocities", data = Vel)
part0.create_dataset("InternalEnergy", data = Uthermal)

## close file
IC.close()

sys.exit(0)
