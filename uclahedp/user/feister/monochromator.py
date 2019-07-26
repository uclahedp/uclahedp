#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
monochromator.py: Routines for analysis of monochromator data

Written by Jessica Pilgrim and Scott Feister, July 2019
"""

import os
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

if __name__ == "__main__":
    h5name = "run006_lens_scan_07-25-2019.hdf5"
    workdir = "/home/scott/myouts/RTFeedback" # Working directory
    h5fn = os.path.join(workdir, h5name)
            
    with h5py.File(h5fn, 'r') as f:
        nsteps = f.attrs["RUN ITERATIONS"]
        
        # # Monochromator PMT
        # data = f['/RESOURCE 133/CHANNEL 0/TRACE']
        # tgv = np.arange(np.size(data,1)) * data.attrs['WAVEFORM DT']
        # shotgv = np.arange(np.size(data,0))

        # xr = tgv[[0,-1]]*1e6
        # yr = shotgv[[0, -1]]
        # fig = plt.figure(1)
        # fig.clear()
        # ax = fig.add_subplot(111)
        # im = ax.pcolorfast(xr, yr, -data[...],
                           # norm=colors.LogNorm(vmin=1e-1, vmax=1e1),
                           # cmap='viridis')
        # cb = fig.colorbar(im, label='-Volts')
        # #ax.set_xlim([8, 12])
        # #ax.set_ylim([0, 300])
        # ax.set_title("Monochromator PMT")
        # ax.set_xlabel("Time (us)")
        # ax.set_ylabel("Shot number")
        
        # Monochromator Timing Diode
        data = f['/RESOURCE 133/CHANNEL 1/TRACE']
        tgv = np.arange(np.size(data,1)) * data.attrs['WAVEFORM DT']
        shotgv = np.arange(np.size(data,0))
        
        # fig = plt.figure(1)
        # fig.clear()

        # xr = tgv[[0,-1]]*1e6
        # yr = shotgv[[0, -1]]
        # fig = plt.figure(1)
        # fig.clear()
        # ax = fig.add_subplot(111)
        # im = ax.pcolorfast(xr, yr, data[...],
                           # norm=colors.Normalize(vmin=-5, vmax=5),
                           # cmap='RdBu')
        # cb = fig.colorbar(im, label='Volts')
        # ax.set_xlim([2.3, 2.6])
        # #ax.set_ylim([0, 300])
        # ax.set_title("Monochromator TDiode")
        # ax.set_xlabel("Time (us)")
        # ax.set_ylabel("Shot number")

        ix0 = np.argmax(np.gradient(data, axis=1), axis=1)
        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(shotgv, ix0)

        print("Max value at " + str(np.argmax(ix0)))
        # print("TGV: " + str(tgv.shape))
        # print("data: " + str(data.shape))
        
        # fig = plt.figure(2)
        # ax = fig.add_subplot(121)
        # ax.plot(tgv*1e6, np.gradient(data[33,:]), '.')
        # ax.set_xlim([2.38,2.395+0.05])

        # ax = fig.add_subplot(122)
        # ax.plot(tgv*1e6, data[33,:], '.')
        # ax.set_xlim([2.38,2.395+0.05])

        print(tgv[1194]*1e6)
        plt.tight_layout()
        fig.savefig(os.path.join(workdir, "MC1.png"), dpi=150)
        print("Analysis complete.")
