""" @package ./examples/current_sheet_2d/create.py
Code that creates 2d current sheet initial conditions

literature reference: Gardiner and Stone (2005), JCoPh.205..509G

created by Rainer Weinberger, last modified 13.03.2019 -- comments welcome
"""

""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/current_sheet_2d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(2.0) # quadratic box
CellsPerDimension = IntType(256)


## parameters
density_0 = FloatType(1.0)
velocity_0 = FloatType(0.0) ## velocity boost
vpert = 0.1
pressure_0 = 0.1
gamma = FloatType(5.0/3.0)
gamma_minus_one = FloatType(gamma - 1.0)
b0 = FloatType(1.0)


""" set up grid """
NumberOfCells = CellsPerDimension*CellsPerDimension
dx = Boxsize / FloatType(CellsPerDimension)
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx
Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
xx, yy = np.meshgrid(Grid1d, Grid1d)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = xx.reshape(NumberOfCells)
Pos[:,1] = yy.reshape(NumberOfCells)


""" set up magnetohydrodynamical quantitites """
## mass insetad of density
Mass = np.full(NumberOfCells, density_0*dx*dx, dtype=FloatType)
## velocity
Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
Velocity[:,0] = velocity_0 + vpert * np.sin(np.pi * Pos[:,1])
Velocity[:,1] = 0.5 * velocity_0
## specific internal energy
Uthermal = np.full(NumberOfCells, (pressure_0/density_0/gamma_minus_one), dtype=FloatType)
## magnetic field strength
Bfield = np.zeros([NumberOfCells, 3], dtype=FloatType)
Bfield[:,1] = b0
i_reverse, = np.where( (Pos[:,0] > 0.5) & (Pos[:,0] < 1.5) )
Bfield[i_reverse, 1] *= -1.0


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
part0.create_dataset("Masses", data = Mass)
part0.create_dataset("Velocities", data = Velocity)
part0.create_dataset("InternalEnergy", data = Uthermal)
part0.create_dataset("MagneticField", data = Bfield)

## close file
IC.close()

exit(0)
