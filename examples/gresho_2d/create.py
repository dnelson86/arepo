""" @package ./examples/Gresho_2d/create.py
Code that creates 2d Gresho vortex initial conditions

created by Rainer Weinberger, last modified 21.02.2019 -- comments welcome
"""

#### load libraries
import sys    # system specific calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/Gresho_2d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(1.0)
CellsPerDimension = IntType(40)

## parameters
density_0 = 1.0
velocity_0 = 3.0 ## bulk velocity
gamma = 5.0/3.0
gamma_minus_one = gamma - 1.0

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
## mass insetad of density
Density = np.full(NumberOfCells, density_0, dtype=FloatType) 
## different zones
i1, = np.where( Radius < 0.2 )
i2, = np.where( (Radius >= 0.2) & (Radius < 0.4) )
i3, = np.where( Radius >= 0.4 )

## velocity
RotationVelocity = np.zeros( NumberOfCells, dtype=FloatType )
RotationVelocity[i1] = 5.0 * Radius[i1]
RotationVelocity[i2] = 2.0 - 5.0 * Radius[i2]
RotationVelocity[i3] = 0.0
Vel = np.zeros([NumberOfCells, 3], dtype=FloatType )
i_all_but_central = np.arange(1,NumberOfCells)
Vel[i_all_but_central,0] = RotationVelocity[i_all_but_central] * yPosFromCenter[i_all_but_central] / Radius[i_all_but_central]
Vel[i_all_but_central,1] = -RotationVelocity[i_all_but_central] * xPosFromCenter[i_all_but_central] / Radius[i_all_but_central]
Vel[:,0] += velocity_0

## specific internal energy
Pressure = np.zeros( NumberOfCells, dtype=FloatType)
Pressure[i1] = 5.0 + 12.5 * Radius[i1] * Radius[i1]
Pressure[i2] = 9.0 + 12.5 * Radius[i2] * Radius[i2] - 20.0 * Radius[i2] + 4.0 * np.log(Radius[i2] / 0.2)
Pressure[i3] = 3.0 + 4.0 * np.log(2.0)
Uthermal = Pressure / density_0 / gamma_minus_one

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
