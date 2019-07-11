""" @package ./examples/galaxy_merger_star_formation_3d/check.py
Code that checks results of 3d galaxy merger problem

created by Rainer Weinberger, last modified 09.03.2019
"""


""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    ## plot stuff
from scipy.interpolate import interp1d    ## inetrpolation
plt.rcParams['text.usetex'] = True

createReferenceSolution = False
makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("galaxy_merger_star_formation_3d: checking simulation output in directory " + simulation_directory) 



FloatType = np.float64  # double precision: np.float64, for single use np.float32
Gcgs = 6.6738e-8


""" compare to reference output """
data = np.loadtxt(simulation_directory+"/output/sfr.txt")
if createReferenceSolution:
    np.savetxt(simulation_directory+"/sfr_reduced_ref.txt", data[::50,:])

data_ref = np.loadtxt(simulation_directory+"/sfr_reduced_ref.txt")
time_ref = data_ref[:,0]
mstar_ref =  data_ref[:,5]

time = data[:,0]
mstar =  data[:,5]

time_probe = np.linspace(0.5,2.8,10)
delta_mass_stars = np.zeros(time_probe.shape)
avg_mass_stars = np.zeros(time_probe.shape)

for i, t in enumerate(time_probe):
    delta_mass_stars[i]  = np.interp(t, time, mstar)
    delta_mass_stars[i] -= np.interp(t, time_ref, mstar_ref)
    avg_mass_stars[i]  = 0.5 * np.interp(t, time, mstar)
    avg_mass_stars[i] += 0.5 * np.interp(t, time_ref, mstar_ref)


if makeplots:
    fig, ax = plt.subplots(2,3, figsize=np.array([6.9,5.35]), sharex=True, sharey=True )
    fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.1,right=0.98, bottom=0.1, top=0.98)

    ax[0, 0].set_xticks([300,320,340])
    ax[0, 0].set_yticks([300,325,350])

for i_plot, i_snap in enumerate([10,11,12,21,22,24]):
    filename = '/output/snap_%03d.hdf5'%i_snap
    print(simulation_directory+filename)
    try:
        data = h5py.File(simulation_directory+filename, "r")
    except:
        break
    ## figure
    pos1 = np.array( data['PartType1']['Coordinates'], dtype=np.float64 )
    pos2 = np.array( data['PartType2']['Coordinates'], dtype=np.float64 )
    pos3 = np.array( data['PartType3']['Coordinates'], dtype=np.float64 )
    pos4 = np.array( data['PartType4']['Coordinates'], dtype=np.float64 )
    
    ix = np.int(i_plot % 3)
    iy = np.int(i_plot / 3)
    
    if makeplots:
        ax[iy,ix].scatter(pos2[:, 0], pos2[:, 1], marker='.', c='k', s=0.05, alpha=0.5, rasterized=True, zorder=2)
        ax[iy,ix].scatter(pos3[:, 0], pos3[:, 1], marker='.', c='k', s=0.05, alpha=0.5, rasterized=True, zorder=2)
        ax[iy,ix].scatter(pos4[:, 0], pos4[:, 1], marker='.', c='red', s=0.05, alpha=0.2, rasterized=True, zorder=3)
        ax[iy,ix].set_aspect('equal')
        ax[iy,ix].set_xlim([295,355])
        ax[iy,ix].set_ylim([290,360])
    

if makeplots:
    ## trajectory
    filename = '/output/snap_%03d.hdf5'%0
    data = h5py.File(simulation_directory+filename, "r")
    pos2 = np.array( data['PartType2']['Coordinates'], dtype=np.float64 )
    ids2 = np.array( data['PartType2']['ParticleIDs'], dtype=np.int32 )
    i_gal1 = np.where(pos2[:,0] < 320)[0]
    i_gal2 = np.where(pos2[:,0] > 330)[0]
    id_gal1 = ids2[i_gal1]
    id_gal2 = ids2[i_gal2]
    
    
    pos_gal1 = np.zeros([32,3], dtype=np.float64)
    pos_gal2 = np.zeros([32,3], dtype=np.float64)
    
    for i_snap in np.arange(32):
        filename = '/output/snap_%03d.hdf5'%i_snap
        try:
            data = h5py.File(simulation_directory+filename, "r")
        except:
            break
        pos2 = np.array( data['PartType2']['Coordinates'], dtype=np.float64 )
        ids2 = np.array( data['PartType2']['ParticleIDs'], dtype=np.int32 )
        ids, i_gal1, dummy = np.intersect1d(ids2, id_gal1, return_indices=True)
        ids, i_gal2, dummy = np.intersect1d(ids2, id_gal2, return_indices=True)
        
        pos_gal1[i_snap,:] = np.average(pos2[i_gal1,:], axis=0)
        pos_gal2[i_snap,:] = np.average(pos2[i_gal2,:], axis=0)
        
    p1x = interp1d(np.arange(32), pos_gal1[:,0], kind='cubic')
    p1y = interp1d(np.arange(32), pos_gal1[:,1], kind='cubic')
    p2x = interp1d(np.arange(32), pos_gal2[:,0], kind='cubic')
    p2y = interp1d(np.arange(32), pos_gal2[:,1], kind='cubic')
    
    
    interp = np.linspace(0,10,1000)
    ax[0,0].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[0,0].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    interp = np.linspace(0,11,1000)
    ax[0,1].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[0,1].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    interp = np.linspace(0,12,1000)
    ax[0,2].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[0,2].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    interp = np.linspace(11,21,1000)
    ax[1,0].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[1,0].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    interp = np.linspace(12,22,1000)
    ax[1,1].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[1,1].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    interp = np.linspace(14,24,1000)
    ax[1,2].plot(p1x(interp), p1y(interp), 'b:', alpha=0.5, lw=4, zorder=1)
    ax[1,2].plot(p2x(interp), p2y(interp), 'b--', alpha=0.5, lw=4, zorder=1)
    
    ax[0,0].set_ylabel("[kpc]")
    ax[1,0].set_ylabel("[kpc]")
    ax[1,0].set_xlabel("[kpc]")
    ax[1,1].set_xlabel("[kpc]")
    ax[1,2].set_xlabel("[kpc]")
    
    if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
    fig.savefig(simulation_directory+'/plots/stars_evolution.pdf')
    
    fig=plt.figure()
    ax = plt.axes([0.15,0.15,0.8,0.8])
    ax.plot(time, mstar * 1e10/0.6774, label="new")
    ax.plot(time_ref, mstar_ref * 1e10/0.6774, label="ref")
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"Stellar mass formed [M$_\odot$]")
    ax.legend(loc=2)
    fig.savefig(simulation_directory+"/plots/Mstar_time.pdf")


""" check if deviations within tolerance """
abs_tolerance = 3.e8 / (1.e10 / 0.6774)
if np.max( np.abs(delta_mass_stars ) ) > abs_tolerance:
    print("Error: stellar mass deviaties from reference more than %g (Msun/Mpc^3): "% abs_tolerance )
    print("times:")
    print(time_probe)
    print("relative deviation:")
    print(delta_mass_stars / avg_mass_stars)
    print("absolute deviation (1e6 Msun):")
    print(delta_mass_stars*1.0e4/0.6774)
    sys.exit(1)

""" if everything is ok """
sys.exit(0)
