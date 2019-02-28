#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:42:19 2019

@author: peter
"""
import numpy as np
import h5py

def spectrogram(f, dt, pad=None, width=None, increment=None, window=None, inTimeUnits=None):
    
    nf = len(f)
    time = np.arange(0, nf)*dt
    
    
    if pad is None:
        pad = True

    
    if inTimeUnits is None:
        inTimeUnits = False
        
        
    if width is None:
        if inTimeUnits:    
            width = dt
        else:
            width = int(nf/100)
    
    if increment is None:
        increment = 0.01*width
        
        
    if window is None:
        window = 'hanning'
        

    if inTimeUnits:
        width_ind = width/dt
        inc_ind = increment/dt
        
        ta = np.argmin( np.abs(  ) )
        
    else:
        width_ind = width
        inc_ind = increment
        
        
    if increment > width:
        raise ValueError("Increment must be smaller than width!")
        
        
        
    if window == 'hanning':
        window = np.hanning(width_ind)
    elif window == 'hamming':
        window = np.hamming(width_ind)
    elif window == 'blackman':
        window = np.blackman(width_ind)
        
    
    
    startind = nf
    stopind = 2*nf-1
    
    plt.plot( f)
    plt.show()
    

    if pad:
        #Pad the array
        temp = np.zeros( (3*nf) )
        temp[0:nf] = np.random.normal( size=(nf) )*np.max(f)*1e-4 + f[0]
        temp[nf:2*nf] = f
        temp[2*nf:3*nf] = np.random.normal( size=(nf) )*np.max(f)*1e-4 + f[-1]
        f = temp
       
        #Pad the time array
        temp = np.zeros( (3*nf) )
        temp[0:nf] = -np.flip(time)
        temp[nf:2*nf] = time
        temp[2*nf:3*nf] = time + time[-1]
        time = temp
        
        nf = len(f)
        
        del(temp)
        
        
    plt.plot(time,f)
    plt.show()
    

    #TODO: Fix this so it incremends in time properly using 'increment'
    window_start = np.arange(nwindows, 
    window_stop = np.arange(startind + int(width_ind/2.0), stopind + int(width_ind/2.0), width_ind)
    window_center = np.arange(startind, stopind, width_ind)
    nwindows = len(window_start)
    window_time = time[window_center]
    
    print('Nwindows: ' + str(nwindows))
    
    nyquist = int(width_ind/2.0)
    freq = np.fft.rfftfreq(width_ind, d=dt)
    nfreq = len(freq)
    
    spectrum = np.zeros( (nwindows, nfreq) )
    
    
    for i in range(nwindows):
    
        #print( str(window_start[i]) + ':' + str(window_stop[i]) )
        f_fft = np.fft.rfft( f[window_start[i]:window_stop[i]]*window )
        spectrum[i,:] =  np.abs(f_fft)**2
        
    return window_time, freq, spectrum

    
    
if __name__ == '__main__':
    f = '/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run40_LAPD7_full.hdf5'
    with h5py.File(f, 'r') as f:
        arr = f['data'][0,:,1]
        time = f['time'][:]
        dt = np.mean(np.gradient(time))
        time, freq, spectrum = spectrogram(arr, dt, width=150, increment=1)
        
        time = time*1e6
        
        fig, ax = plt.subplots( figsize = [6, 6])
        cplot = ax.imshow(pow(spectrum, 0.15), cmap='bone',
                          extent = [time[0], time[-1],freq[0], freq[-1]],
                          aspect='auto',origin='lower')
        ax.set_ylim( (0, 1e7) )
        ax.set_xlim( (0, 30) )
