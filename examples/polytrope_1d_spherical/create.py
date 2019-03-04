""" @package examples/polytrope_1d_spherical/create.py
Code that creates 1d polytrope test problem;
supposed to be as siple as possible and self-contained

created by Rainer Weinberger, last modified 04.03.2019
"""


""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/polytrope_1d_spherical/create.py: creating ICs in directory " + simulation_directory)


""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(1.0)
NumberOfCells = IntType(256)

## initial state
G = FloatType(1.0)
density_0 = FloatType(1.0)
PolytropicIndex = FloatType(1.0)
velocity_0 = FloatType(0.0)
gamma = FloatType( (1.0 + PolytropicIndex) / PolytropicIndex)  ## note: this has to be consistent with the parameter settings for Arepo
pressure_0 = FloatType(1.)/gamma


""" set up grid: uniform radial 1d grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)
## position of first and last cell
pos_first, pos_last = FloatType(0.1) + FloatType(0.5) * dx, Boxsize - FloatType(0.5) * dx

## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = FloatType(4.) * np.pi / FloatType(3.) * ( (Pos[:,0]+FloatType(0.5)*dx)**3 - (Pos[:,0]-FloatType(0.5)*dx)**3 )


""" set up hydrodynamical quantitites """
Rscale = FloatType(0.8) * Boxsize / np.pi
Radius = Pos[:,0]
theta = np.sinc(Radius/Rscale/np.pi) ## np.sinc(x) = np.sin(np.pi x)/(np.pi x) or 1 if x=0
theta[theta<1.e-12] = FloatType(1.e-12)
K = Rscale**2 * FloatType(4.0) * np.pi * G / ( (PolytropicIndex+FloatType(1.0) ) * density_0**(FloatType(1.0)/PolytropicIndex - FloatType(1.0) ) )
## explicit n=1 case
K = Rscale**2 * FloatType(2.0) * np.pi * G

## hydrodynamic quantities
Density = density_0 * theta**PolytropicIndex
Velocity = np.zeros([NumberOfCells, 3], dtype=FloatType)
Velocity[:,0] = velocity_0
Pressure = K * Density**gamma
Uthermal = Pressure / Density / (gamma - FloatType(1.0) )

## density in mass filed in input
Mass = Density


""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/create file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
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

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)

## close file
IC.close()


""" normal exit """
sys.exit(0) 
