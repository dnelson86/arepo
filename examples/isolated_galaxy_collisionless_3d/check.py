""" @package ./examples/isolated_galaxy_collisionless_3d/check.py
Code that checks results of 3d isolated collisionless galaxy problem

created by Rainer Weinberger, last modified 10.03.2019
"""


""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    ## plot stuff
from scipy.interpolate import interp1d    ## inetrpolation
plt.rcParams['text.usetex'] = True
from matplotlib.colors import LogNorm

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("examples/isolated_galaxy_collisionless_3d/check.py: checking simulation output in directory " + simulation_directory) 


FloatType = np.float64  # double precision: np.float64, for single use np.float32
Gcgs = 6.6738e-8
status=0


rGrid = np.linspace(1,34,100)
mEnc1Grid_ref = np.zeros(rGrid.shape)
mEnc2Grid_ref = np.zeros(rGrid.shape)
mEnc3Grid_ref = np.zeros(rGrid.shape)


""" loop over all output files """
i_file = 0
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    
    """ get simulation data """
    time = FloatType(data["Header"].attrs["Time"])
    Pos1 = np.array(data["PartType1"]["Coordinates"], dtype = FloatType)
    Mass1 = np.array(data["PartType1"]["Masses"], dtype=FloatType)
    Pos2 = np.array(data["PartType2"]["Coordinates"], dtype = FloatType)
    Mass2 = np.array(data["PartType2"]["Masses"], dtype=FloatType)
    Pos3 = np.array(data["PartType3"]["Coordinates"], dtype = FloatType)
    Mass3 = np.array(data["PartType3"]["Masses"], dtype=FloatType)
    
    
    Center = np.zeros(3) ## center defined as center of disk
    Center[0] = np.average(Pos2[:,0])
    Center[1] = np.average(Pos2[:,1])
    Center[2] = np.average(Pos2[:,2])
    
    print('center: %g %g %g'%(Center[0], Center[1], Center[2]))
    
    xPosFromCenter1 = (Pos1[:,0] - Center[0])
    yPosFromCenter1 = (Pos1[:,1] - Center[1])
    zPosFromCenter1 = (Pos1[:,2] - Center[2])
    
    
    xPosFromCenter2 = (Pos2[:,0] - Center[0])
    yPosFromCenter2 = (Pos2[:,1] - Center[1])
    zPosFromCenter2 = (Pos2[:,2] - Center[2])
    
    
    xPosFromCenter3 = (Pos3[:,0] - Center[0])
    yPosFromCenter3 = (Pos3[:,1] - Center[1])
    zPosFromCenter3 = (Pos3[:,2] - Center[2])
    
    
    Radius1 = np.sqrt( xPosFromCenter1**2 + yPosFromCenter1**2 + zPosFromCenter1**2 )
    Radius2 = np.sqrt( xPosFromCenter2**2 + yPosFromCenter2**2 + zPosFromCenter2**2 )
    Radius3 = np.sqrt( xPosFromCenter3**2 + yPosFromCenter3**2 + zPosFromCenter3**2 )
    
    
    """ disk height """
    h_disp = np.std(zPosFromCenter2)
    if h_disp > 1.2:
        status = 1
        print("ERROR: %d disk height: %g"% (i_file, h_disp) )
    
    
    """ rotation curve / enclosed mass """
    ## sort by radius
    i_sort1 = np.argsort(Radius1)
    i_sort2 = np.argsort(Radius2)
    i_sort3 = np.argsort(Radius3)
    
    
    ## calculate enclosed mass
    mEnc1 = np.cumsum(Mass1[i_sort1])
    mEnc2 = np.cumsum(Mass2[i_sort2])
    mEnc3 = np.cumsum(Mass3[i_sort3])
    
    
    #interpolate on uniform grid
    mEnc1Grid = np.interp( rGrid, Radius1[i_sort1], mEnc1 )
    mEnc2Grid = np.interp( rGrid, Radius2[i_sort2], mEnc2 )
    mEnc3Grid = np.interp( rGrid, Radius3[i_sort3], mEnc3 )
    mEncTotGrid = mEnc1Grid + mEnc2Grid + mEnc3Grid
    
    
    ## calculate rotation curve
    vCirc1Grid = np.sqrt(mEnc1Grid * Gcgs / rGrid * 6.4459e11)
    vCirc2Grid = np.sqrt(mEnc2Grid * Gcgs / rGrid * 6.4459e11)
    vCirc3Grid = np.sqrt(mEnc3Grid * Gcgs / rGrid * 6.4459e11)
    vTotGrid =  np.sqrt(mEncTotGrid * Gcgs / rGrid * 6.4459e11)
    
    
    if makeplots:
        ## figure of positions
        fig = plt.figure(figsize = np.array([6.9,9.2]))
        ax = []
        ax.append(plt.axes([0.1,0.41,0.75,0.56]) )
        ax.append(plt.axes([0.1,0.07,0.75,0.282]) )
        cax = plt.axes([0.86,0.07,0.04,0.90])
        
        xext = 15.0
        
        binsx = np.linspace(-xext,xext,128)
        binsy = np.linspace(-xext,xext,128)
        binsz = np.linspace(-0.5*xext,0.5*xext,64)
        
        H1, xgrid, ygrid = np.histogram2d(xPosFromCenter2, yPosFromCenter2, bins=[binsx, binsy], density=True)
        H2, xgrid, ygrid = np.histogram2d(xPosFromCenter3, yPosFromCenter3, bins=[binsx, binsy], density=True)
        He1, xgrid, ygrid = np.histogram2d(zPosFromCenter2, yPosFromCenter2, bins=[binsz, binsy], density=True)
        He2, xgrid, ygrid = np.histogram2d(zPosFromCenter3, yPosFromCenter3, bins=[binsz, binsy], density=True)
        
        ext=(-xext,xext,-xext,xext)
        ax[0].imshow(H1+H2, \
                    norm=LogNorm(vmin=3e-4,vmax=3e-1), \
                    cmap='Greys', \
                    origin='lower', \
                    extent=ext \
                    )
        
        ext=(-xext,xext,-0.5*xext,0.5*xext)
        img = ax[1].imshow(He1+He2, \
                    norm=LogNorm(vmin=3e-4,vmax=3e-1), \
                    cmap='Greys', \
                    origin='lower', \
                    extent=ext \
                    )
        
        ax[0].axis('equal')
        ax[1].axis('equal')
        ax[0].set_xticks([-10,-5,0,5,10])
        ax[0].set_yticks([-10,-5,0,5,10])
        ax[1].set_xticks([-10,-5,0,5,10])
        ax[1].set_yticks([-5,0,5])
        
        ax[0].set_xlabel('[kpc]')
        ax[0].set_ylabel('[kpc]')
        ax[1].set_xlabel('[kpc]')
        ax[1].set_ylabel('[kpc]')
        
        fig.colorbar(img, cax=cax)
        
        if not os.path.exists( simulation_directory+"/plots" ):
            os.mkdir( simulation_directory+"/plots" )
        fig.savefig(simulation_directory+"/plots/positions_%03d.pdf"%i_file)
        
        ## figure of rotation curve
        fig = plt.figure()
        ax = plt.axes([0.1, 0.1, 0.85, 0.85])
        ax.plot( rGrid, vCirc1Grid )
        ax.plot( rGrid, vCirc2Grid )
        ax.plot( rGrid, vCirc3Grid )
        ax.plot( rGrid, vTotGrid )
        ax.set_xlim(0,35)
        ax.set_ylim(0,300)
        ax.set_xlabel(r"radius")
        ax.set_ylabel(r"v$_c$")
        fig.savefig(simulation_directory+"/plots/rotation_curve_%03d.pdf"%i_file)
    
    
    ## comparision to first snapshot (ICs)
    if i_file == 0:
        mEnc1Grid_ref[:] = mEnc1Grid[:]
        mEnc2Grid_ref[:] = mEnc2Grid[:]
        mEnc3Grid_ref[:] = mEnc3Grid[:]
    else:
        tolerance = 1.0
        delta_mEnc = mEnc1Grid_ref - mEnc1Grid
        if np.max( abs(delta_mEnc) ) > tolerance:
            status = 1
            print("ERROR: i_file: %d: enclused mass difference in halo (part1) exceeds tolerance %g" % (i_file, tolerance) )
            print(delta_mEnc)
        
        
        delta_mEnc = mEnc2Grid_ref - mEnc2Grid
        if np.max( abs(delta_mEnc) ) > tolerance:
            status = 1
            print("ERROR: i_file: %d: enclused mass difference in disc (part2) exceeds tolerance %g" % (i_file, tolerance) )
            print(delta_mEnc)
        
        
        delta_mEnc = mEnc3Grid_ref - mEnc3Grid
        if np.max( abs(delta_mEnc) ) > tolerance:
            status = 1
            print("ERROR: i_file: %d: enclused mass difference in bulge (part3) exceeds tolerance %g" % (i_file, tolerance) )
            print(delta_mEnc)
    i_file += 1


## if everything is ok: 0; else: 1
sys.exit(status)
