#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rtfeedback.py: Real-time feedback for July 2019 Run

Analyzes HDF5 file and plots data. Wrap in a bash sheel while loop to continuously attempts to update file and replot.

Created by Scott Feister on Wed Jul 24 14:44:13 2019
"""

import os
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

if __name__ == "__main__":
    print("Beginning synchonization.")
    h5name = "run007_lens_scan_07-25-2019.hdf5"
    workdir = "/home/scott/myouts/RTFeedback" # Working directory
    h5fn = os.path.join(workdir, h5name)
    
    #os.system('rsync -v sfeister@128.97.13.200:"/home/shared/Phoenix\ Terminal/Data/HRR/' + h5name + '" "' + workdir + '/"')
    print("Beginning analysis.")
        
    with h5py.File(h5fn, 'r') as f:
        nsteps = f.attrs["RUN ITERATIONS"]
        
        # Monochromator
        data = f['/RESOURCE 133/CHANNEL 0/TRACE']
        tgv = np.arange(np.size(data,1)) * data.attrs['WAVEFORM DT']
        shotgv = np.arange(np.size(data,0))

        xr = tgv[[0,-1]]*1e6
        yr = shotgv[[0, -1]]
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(221)
        im = ax.pcolorfast(xr, yr, -data[...],
                           norm=colors.LogNorm(vmin=1e-1, vmax=1e1),
                           cmap='viridis')
        cb = fig.colorbar(im, label='-Volts')
        #ax.set_xlim([8, 12])
        #ax.set_ylim([0, 300])
        ax.set_title("Monochromator PMT")
        ax.set_xlabel("Time (us)")
        ax.set_ylabel("Shot number")

        # Magnetic field components
        titles = ['Bx', 'By', 'Bz'] # B components
        vmaxes = np.array([0.5,0.5,0.5])*1e3 # For each B component, colorbar +/- limits in millivolts
        chans = [1,2,3] # For each B component, colorbar limit
        for i in range(3):
            data = f['/RESOURCE 123/CHANNEL ' + str(chans[i]) + '/TRACE']
            tgv = np.arange(np.size(data,1)) * data.attrs['WAVEFORM DT']
            shotgv = np.arange(np.size(data,0))

            xr = tgv[[0,-1]]*1e6
            yr = shotgv[[0, -1]]
            
            ax = fig.add_subplot(2, 2, i + 2)
            im = ax.pcolorfast(xr, yr, 1e3*data[...],
                               norm=colors.Normalize(vmin=-vmaxes[i], vmax=vmaxes[i]),
                               cmap='RdBu')
            cb = fig.colorbar(im, label='mV')
            ax.set_title(titles[i])
            ax.set_xlabel("Time (us)")
            ax.set_ylabel("Shot number")
            
        plt.tight_layout()
        fig.savefig(os.path.join(workdir, "RT1.png"), dpi=150)
        print("Analysis complete.")
            
            
