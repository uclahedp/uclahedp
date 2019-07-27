#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
monochromator.py: Quick plot of monochromator data.

Analyzes the timing diode and subtracts that from the overall traces.
Reads the specified HDF5 RAW file, located in the workdir.
Generates plots including the velocity distribution.
Plots are saved into the workdir as 'RAW.png', 'VEL.png', 'POS.png'.

Written by Scott Feister on July 25, 2019
"""

import os
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def indexed_roll(a, shift):
    """ Rolls the 2D array a along axis 1
    with a variable roll quantity based on the 1D array shift
    """
    if a.shape[0] != len(shift):
        raise Exception("Dim0 of a does not match length of shift array")

    b = np.zeros(a.shape)

    for i in range(len(shift)):
        b[i,:] = np.roll(a[i,:], shift[i])
        
        # Replace areas with NAN (this doesn't seem to work)
        if shift[i] < 0:
            b[i,b.shape[1]-shift[i]:] = np.nan
        if shift[i] > 0:
            b[i,:shift[i]] = np.nan            
    
    return b

if __name__ == "__main__":
    workdir = "/home/scott/myouts/RTFeedback" # Working directory; both where the HDF5 file is located and where the plot is saved.
    h5name = "run006_lens_scan_07-25-2019.hdf5" # HDF5 file to analyze (located in workdir)
    h5fn = os.path.join(workdir, h5name)
            
    with h5py.File(h5fn, 'r') as f:     
        portTarg = 13 # Target port number at LAPD
        portPMT = 14 # PMT port number at LAPD
        dist2targ = np.abs(portPMT - portTarg) * 0.325 # Distance from target to PMT, in meters
        print("Distance from PMT port to target port: " + str(dist2targ*1e3) + " mm")
        
        ## Extract the lens positions array
        posgv = f['/RESOURCE 145/CHANNEL 0/POSITION'][...] # Lens positions, in mm

        ## Monochromator Timing Diode: Find time-zero T0 shift (for each shot)
        # Define t0 as point of maximum time derivative
        tdata = f['/RESOURCE 133/CHANNEL 1/TRACE'] # Timing diode
        ix0 = np.argmax(np.gradient(tdata, axis=1), axis=1) # 1D array, number of indices to shift each trace
        
        ## Monochromator PMT
        data = f['/RESOURCE 133/CHANNEL 0/TRACE']
        tgv = np.arange(np.size(data,1)) * data.attrs['WAVEFORM DT'] # Time axis values in seconds
        shotgv = np.arange(np.size(data,0)) # Shotnumber axis values
        
        A = -indexed_roll(data[...],-ix0) # Shift data indices to adjust to timing diode
        
        ## MAKE PLOTS
        # Position graph
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(111)
        ax.plot(posgv)
        ax.set_ylabel("Lens position (mm)")
        ax.set_xlabel("Step number")
        plt.grid(b=True,which='both')
        fig.savefig(os.path.join(workdir, "POS.png"), dpi=150)

        # PMT Raw Graph
        xr = tgv[[0,-1]]*1e6 # [xmin, xmax]
        yr = shotgv[[0, -1]] # [ymin, ymax]
        
        fig = plt.figure(2)
        fig.clear()
        ax = fig.add_subplot(221)
        im = ax.pcolorfast(xr, yr, tdata[...], vmin=0, cmap='plasma')
        cb = fig.colorbar(im, label='Volts')
        ax.set_title("Raw Timing Diode")
        ax.set_xlabel("Time (us)")
        ax.set_ylabel("Shot number")

        ax = fig.add_subplot(222)
        ax.plot(shotgv, ix0)
        ax.set_title("Timing Diode Offsets")
        ax.set_ylabel("Indices offset")
        ax.set_xlabel("Shot number")
        
        ax = fig.add_subplot(223)
        im = ax.pcolorfast(xr, yr, -data[...],
                           norm=colors.LogNorm(vmin=1e-1, vmax=1e1),
                           cmap='viridis')
        cb = fig.colorbar(im, label='-Volts')
        ax.set_title("Raw PMT")
        ax.set_xlabel("Time (us)")
        ax.set_ylabel("Shot number")

        ax = fig.add_subplot(224)
        im = ax.pcolorfast(xr, yr, A, #-data[...],
                           norm=colors.LogNorm(vmin=1e-1, vmax=1e1),
                           cmap='viridis')
        cb = fig.colorbar(im, label='-Volts')
        ax.set_title("Time-shifted PMT")
        ax.set_xlabel("Time (us)")
        ax.set_ylabel("Shot number")

        plt.tight_layout()
        fig.savefig(os.path.join(workdir, "RAW.png"), dpi=150)
        
        # Velocity Graph
        smax = np.max(ix0) # Maximum index shift (due to timing diode)
        velgv = dist2targ / tgv[1:-smax][::-1] # Velocity axis in m/s, shortened and reversed
        A2 = A[:,1:-smax][:,::-1] * velgv # Amplitudes, shortened and reversed, then multiplied by velocity (to account transforming of axis; noting that dx = v dt)
        
        fig = plt.figure(3)
        fig.clear()
        ax = fig.add_subplot(111)
        im = ax.pcolorfast(velgv*1e-3, shotgv, A2[:-1,:-1]*1e-5, cmap='viridis', vmin=0, vmax=7) # The 1e5 is an arbitrary reduction to make smaller A.U. values
        cb = fig.colorbar(im, label='photons/dv (A.U.)')
        ax.set_xlim([0, 400])
        ax.set_title("Velocity Distribution")
        ax.set_xlabel("Velocity (km/s)")
        ax.set_ylabel("Shot number")
        
        fig.savefig(os.path.join(workdir, "VEL.png"), dpi=150)

        print("Analysis complete.")
