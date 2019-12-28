""" Common utilities for AREPO test suite. """
import subprocess
import numpy as np
import glob
import h5py
import os
import os.path


def fail(msg):
    print(' - \033[91m%s\033[0m' % msg)


def success(msg):
    print(' - \033[92m%s\033[0m' % msg)


parameterFileVisOptions = """
PicXpixels 256
PicYpixels 256
PicXaxis 0
PicYaxis 1
PicZaxis 2
PicXmin 0.0
PicXmax MAX_X
PicYmin 0.0
PicYmax MAX_Y
PicZmin 0.0
PicZmax 0.0
TimeBetweenImages 0.05"""


def write_ic_file(fileName, partTypes, boxSize, massTable=None,
                  headerAttrs={}):
    """ Helper to write a HDF5 IC file. partTypes is a dictionary with keys of the form 
    PartTypeX, each of which is its own dictionary of particle fields and ndarrays. 
    boxSize is a scalar float, and massTable a 6-element float array, if specified. 
    If headerAttrs is not an empty, then a dict of header attributes to set. """
    import h5py

    nPartTypes = 6

    with h5py.File(fileName, 'w') as f:
        # write each PartTypeX group and datasets
        for ptName in partTypes.keys():
            g = f.create_group(ptName)
            for field in partTypes[ptName]:
                g[field] = partTypes[ptName][field]

        # set particle counts
        NumPart = np.zeros(nPartTypes, dtype='int32')
        for ptName in partTypes.keys():
            ptNum = int(ptName[-1])
            NumPart[ptNum] = partTypes[ptName]['ParticleIDs'].size

        # create standard header
        h = f.create_group('Header')
        h.attrs['BoxSize'] = boxSize
        h.attrs['NumFilesPerSnapshot'] = 1
        h.attrs['NumPart_ThisFile'] = NumPart
        h.attrs['NumPart_Total'] = NumPart
        h.attrs['NumPart_Total_HighWord'] = np.zeros(nPartTypes, dtype='int32')

        for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[k] = 0.0
        for k in [
                'Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback',
                'DoublePrecision'
        ]:
            h.attrs['Flag_%s' % k] = 0
        for k in headerAttrs:
            h.attrs[k] = headerAttrs[k]

        if massTable is not None:
            h.attrs['MassTable'] = massTable
        else:
            h.attrs['MassTable'] = np.zeros(nPartTypes, dtype='float64')


def visualize_result_2d(path, Lx, Ly):
    """ Helper function to load density_field_NNN projection files and plot a series of PNG frames. """
    import struct
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    vMM = None
    figsize = (16, 14)

    # loop over density projections
    snapPath = os.path.join(path, 'output')
    savePath = os.path.join(path, 'vis')

    nSnaps = len(glob.glob(os.path.join(snapPath, 'density_field_*')))

    for i in range(nSnaps):
        # load
        with open(os.path.join(snapPath, 'density_field_%03d' % i),
                  mode='rb') as f:
            data = f.read()

        # unpack
        nPixelsX = struct.unpack('i', data[0:4])[0]
        nPixelsY = struct.unpack('i', data[4:8])[0]

        nGridFloats = (len(data) - 8) / 4
        grid = struct.unpack('f' * nGridFloats, data[8:])
        grid = np.array(grid).reshape((nPixelsX, nPixelsY))

        if vMM is None: vMM = [grid.min(), grid.max()]  # set on first snap

        # plot
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('%3d' % i)

        plt.imshow(grid.T, extent=[0, Lx, 0, Ly], cmap='magma', aspect=1.0)
        ax.autoscale(False)
        plt.clim(vMM)

        # colorbar
        cax = make_axes_locatable(ax).append_axes('right', size='4%', pad=0.2)
        cb = plt.colorbar(cax=cax)
        cb.ax.set_ylabel('Density')

        fig.tight_layout()
        fig.savefig(os.path.join(savePath, 'density_%d.png' % i))
        plt.close(fig)

    # if ffmpeg exists, make a movie
    cmd = [
        'ffmpeg', '-framerate', '10', '-i', 'density_%d.png', '-pix_fmt',
        'yuv420p', 'movie.mp4', '-y'
    ]
    try:
        output = subprocess.check_output(cmd,
                                         stderr=subprocess.STDOUT,
                                         cwd=os.path.join(path, 'vis'))
    except:
        pass


def visualize_result_1d(path, Lx, radial=False, showICs=True):
    """ Helper function to load and plot a series of snapshots for a simple 1D run. 
    If radial == True, plot 1D profiles as a function of the radial distance from the (square) box center.
    If showICs == True, overplot the (constant) initial condition state on each frame. """
    import matplotlib.pyplot as plt

    figsize = (20, 12)
    snapPath = os.path.join(path, 'plots')
    savePath = os.path.join(path, 'vis')
    icsPath = os.path.join(path, 'ics.hdf5')
    ptName = 'PartType0'

    if not os.path.exists(snapPath):
        os.mkdir(snapPath)

    nSnaps = len(glob.glob(os.path.join(snapPath, 'snap_*.hdf5')))

    if showICs:
        ics = {}
        with h5py.File(icsPath, 'r') as f:
            for k in f[ptName].keys():
                ics[k] = f[ptName][k][()]

        ics['x'] = ics['Coordinates'][:, 0]
        ics['Density'] = ics['Masses']  # fixed assumption if showICs == True
        ics['Pressure'] = ics[
            'Density']**1.4  # fixed gamma assumption if showICs == True

        if radial:
            ics['x'] = np.sqrt((ics['Coordinates'][:, 0] - Lx / 2)**2 +
                               (ics['Coordinates'][:, 1] - Lx / 2)**2)

    for i in range(nSnaps):
        # load
        data = {}
        with h5py.File(os.path.join(snapPath, 'snap_%03d.hdf5' % i), 'r') as f:
            for k in f[ptName].keys():
                data[k] = f[ptName][k][()]
            data['Time'] = f['Header'].attrs['Time']

        # decide x coordinate
        xlabel = 'x' if not radial else 'radius'
        x = data['Coordinates'][:, 0]
        if radial:
            x = np.sqrt((data['Coordinates'][:, 0] - Lx / 2)**2 +
                        (data['Coordinates'][:, 1] - Lx / 2)**2)

        # plot
        fig = plt.figure(figsize=figsize)

        ax = fig.add_subplot(2, 2, 1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Density')
        ax.set_ylim([0.0, 1.1])
        ax.set_title('t = %.4g' % data['Time'])
        ax.plot(x, data['Density'], 'o')
        if showICs:
            ax.plot(ics['x'],
                    ics['Density'],
                    'o',
                    markersize=1.0,
                    color='black')

        ax = fig.add_subplot(2, 2, 2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Pressure')
        ax.set_ylim([0.0, 1.1])
        ax.set_title('t = %.4g' % data['Time'])
        ax.plot(x, data['Pressure'], 'o')
        if showICs:
            ax.plot(ics['x'],
                    ics['Pressure'],
                    'o',
                    markersize=1.0,
                    color='black')

        ax = fig.add_subplot(2, 2, 3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Velocity [x]')
        ax.set_ylim([-0.1, 1.5])
        ax.set_title('t = %.4g' % data['Time'])
        ax.plot(x, data['Velocities'][:, 0], 'o')
        if showICs:
            ax.plot(ics['x'],
                    ics['Velocities'][:, 0],
                    'o',
                    markersize=1.0,
                    color='black')

        ax = fig.add_subplot(2, 2, 4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'Entropy [P/$\rho^\gamma$]')
        ax.set_ylim([0.9, 1.4])
        ax.set_title(r't = %.4g' % data['Time'])
        ax.plot(x, data['Pressure'] / data['Density']**1.4, 'o')
        if showICs:
            ax.plot(ics['x'],
                    ics['Pressure'] / ics['Density']**1.4,
                    'o',
                    markersize=1.0,
                    color='black')

        fig.tight_layout()
        fig.savefig(os.path.join(savePath, 'frame_%d.png' % i))
        plt.close(fig)

    # if ffmpeg exists, make a movie
    cmd = [
        'ffmpeg', '-framerate', '10', '-i', 'frame_%d.png', '-pix_fmt',
        'yuv420p', 'movie.mp4', '-y'
    ]
    try:
        output = subprocess.check_output(cmd,
                                         stderr=subprocess.STDOUT,
                                         cwd=savePath)
    except:
        pass
