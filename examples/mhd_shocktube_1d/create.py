""" @package ./examples/mhd_shocktube_1d/create.py
Code that creates 1d mhd shocktube initial conditions


created by Rainer Weinberger, last modified 12.03.2019 -- comments welcome
"""

""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/mhd_shocktube_1d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32
Boxsize = FloatType(2.5) # quadratic box
NumberOfCells = IntType(400)

alpha = np.pi

## parameters
density_L = FloatType(1.0)
velocity_L = FloatType(0.0)
pressure_L = FloatType(1.0)
b_L = np.array([1.0, 1.0, 0.0], dtype=FloatType)

density_R = FloatType(0.2)
velocity_R = FloatType(0.0)
pressure_R = FloatType(0.2)
b_R = np.array([1.0, np.cos(alpha), np.sin(alpha)], dtype=FloatType)

gamma = FloatType(5.0/3.0)
gamma_minus_one = FloatType(gamma - 1.0)


""" set up grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx

## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = np.full(NumberOfCells, dx, dtype=FloatType)


""" set up magnetohydrodynamical quantitites """
## left state
Mass = np.full(NumberOfCells, density_L*dx, dtype=FloatType)
Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
Velocity[:,0] = velocity_L
Uthermal = np.full(NumberOfCells, (pressure_L/density_L/gamma_minus_one), dtype=FloatType)
Bfield = np.zeros([NumberOfCells, 3], dtype=FloatType)
for dim in np.arange(3):
    Bfield[:,dim] = b_L[dim]

## right state
i_right, = np.where(Pos[:,0] > 1.0)
Mass[i_right] = density_R*dx
Velocity[i_right,0] = velocity_R
Uthermal[i_right] = pressure_R/density_R/gamma_minus_one
for dim in np.arange(3):
    Bfield[:,dim] = b_R[dim]


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
