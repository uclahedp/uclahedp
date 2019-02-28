#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:42:19 2019

@author: peter
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from uclahedp.tests.synthdata import twoWavePackets as twoWavePackets

def spectrogram(f, dt, width=None, increment=None,
                window=None, inTimeUnits=None):

    nf = len(f)

    # Set default values
    if inTimeUnits is None:
        inTimeUnits = False

    if width is None:
        if inTimeUnits:
            width = dt*nf/50
        else:
            width = nf/50

    if increment is None:
        increment = 0.01*width

    if window is None:
        window = 'hanning'

    # If units are given in time units, convert to indices using dt
    if inTimeUnits:
        width_ind = int(np.ceil(width/dt))
        inc_ind = int(increment/dt)

    else:
        width_ind = int(np.ceil(width))
        inc_ind = int(increment)

    # Width_ind must be even (based on way windows are chosen)
    if width_ind % 2 == 1:
        width_ind += 1

    if increment > width:
        raise ValueError("Increment must be smaller than width!")

    if window == 'hanning':
        window = np.hanning(width_ind)
    elif window == 'hamming':
        window = np.hamming(width_ind)
    elif window == 'blackman':
        window = np.blackman(width_ind)
    elif window == 'flattop':
        window = np.ones(width_ind)

    # reate a time vector
    time = np.arange(0, nf)*dt

    # Pad the array
    temp = np.zeros((3*nf))
    temp[0:nf] = np.random.normal(size=(nf))*np.max(f)*1e-4 + f[0]
    temp[nf:2*nf] = f
    temp[2*nf:3*nf] = np.random.normal(size=(nf))*np.max(f)*1e-4 + f[-1]
    f = temp

    # Pad the time array
    temp = np.zeros((3*nf))
    temp[0:nf] = -np.flip(time)
    temp[nf:2*nf] = time
    temp[2*nf:3*nf] = time + time[-1]
    time = temp
    del(temp)
    
    # Compute the center of each window, and the end point indices
    window_center = np.arange(nf, 2*nf-1, inc_ind)
    window_start = window_center - int(width_ind/2.0)
    window_stop = window_center + int(width_ind/2.0)
    

    # Compute the number of windows
    nwindows = len(window_start)
    # Calculate the time cooresponding to each window's center
    window_time = time[window_center]

    # Create a frequency vector
    freq = np.fft.rfftfreq(width_ind, d=dt)
    nfreq = len(freq)

    # Create ouptut array
    spectrum = np.zeros((nwindows, nfreq))
    

    # Loop through each window and compute the FFT
    for i in range(nwindows):
        #print(str(window_start[i]) + '->' + str(window_stop[i]) + ', cent: ' + str(window_center[i]) )
        f_fft = np.fft.rfft(f[window_start[i]:window_stop[i]]*window)
        spectrum[i, :] = np.abs(f_fft)**2

        
    return window_time, freq, spectrum

def _test_spectrogram():
    
    fig, ax = plt.subplots(nrows=2, figsize = [3, 6])
    
    dt, t, arr = twoWavePackets()
    ax[0].plot(t*1e6, arr)
    
    time, freq, spectrum = spectrogram(arr, dt)
    
    ax[1].imshow(pow(spectrum.T, 0.15), cmap='jet',
                          extent = [time[0]*1e6, time[-1]*1e6,freq[0], freq[-1]],
                          aspect='auto',origin='lower')
    ax[1].set_ylim((0, 1e7))

    
    

if __name__ == '__main__':
    #_test_spectrogram()
    
    
    
    f = '/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run40_LAPD7_full.hdf5'
    f = '/Volumes/PVH_DATA/LAPD_Jan2019/FULL/run34_LAPD10_full.hdf5'
    
    zrange = (None,4)
    
    with h5py.File(f) as f:
        arr = f['data'][0,:,1]
        
        print(arr.shape)
     
        t = f['time'][:]
        dt = np.mean(np.gradient(t))
        
        
        time, freq, spectrum = spectrogram(arr, dt, width=1e-5, inTimeUnits=True)
        fig, ax = plt.subplots(nrows=2, figsize = [3, 6])
        ax[0].plot(t*1e6, arr)
        ax[0].set_xlim((0,40))
        cplot = ax[1].imshow(pow(spectrum.T, 0.1), cmap='jet',
                  vmin = zrange[0], vmax=zrange[1],
                  extent = [time[0]*1e6, time[-1]*1e6,freq[0], freq[-1]],
                  aspect='auto',origin='lower')
        ax[1].set_ylim((0, 6e6))
        ax[1].set_xlim((0,40))
        ax[1].ticklabel_format(scilimits=(-3, 3))
        
        cb = fig.colorbar(cplot, orientation='vertical', format = '%.2f')
        
       
        
   
