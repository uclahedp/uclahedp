#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:13:31 2019

@author: peter

***
This program requires ffmpeg be installed in anaconda
To install, run the following command in your anaconda terminal
conda install -c conda-forge ffmpeg
***
"""

import h5py
import numpy as np

import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter


def pickTimeInds(pick_time, timearr):
    output = np.zeros(pick_time.shape, dtype=np.int32)
    for i, t in enumerate(pick_time):
        ti = np.argmin(np.abs(timearr - t))
        output[i] = ti
    return output


def simpleMovie(inds, x, y, arr, savefile, contour=False, cmap=None):
    
    dpi =600
 
    if cmap is None:
        cmap = 'jet'
        
    nframes = len(inds)

    extent = [x[0], x[-1], y[0], y[-1]]

    data = arr[inds, :, :]
    vmin = np.min(data)
    vmax = np.max(data)

    
    metadata = dict(title='Movie test', artist="uclahedp", comment='comment')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    

    fig = plt.figure(figsize=(6,6))
    
    x = np.squeeze(arr[inds[0],:,:])
    cplot = plt.imshow(x, cmap=cmap, aspect='auto',
                       vmin=vmin, vmax=vmax,
                       origin='lower', extent=extent)
    

    with writer.saving(fig, savefile, dpi):
        for i in range(nframes):
            print("Frame: " + str(i))
            x = np.squeeze(arr[inds[i],:,:])
            cplot.set_data(x)
            writer.grab_frame()
            
            
            

if __name__ == '__main__':
    f = '/Volumes/PVH_DATA/LAPD_Jan2019/FULL/run18_LAPD_C6_full.hdf5'
    savefile = '/Volumes/PVH_DATA/LAPD_Jan2019/FULL/simpleMovie.mp4'
    
    with h5py.File(f, 'r') as f:
        x = f['xaxis'][:]
        y = f['yaxis']
        timearr = f['time'][:]
        arr = f['data'][:,:,:,0,0,2]
        
        times = np.linspace(0, 2, num=60)
        tinds = pickTimeInds(times, timearr*1e6)
        
        simpleMovie(tinds, x, y, arr, savefile)
        
        