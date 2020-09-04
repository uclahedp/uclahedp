#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 13:08:55 2019

@author: peter
This file contains filtering routines
"""
from uclahedp.tests import synthdata

import matplotlib.pyplot as plt
import numpy as np

import scipy.signal as signal



def fftFilter(f, dt, band=(None,None), mode='pass', plots=False, plotband=None):
    """
    band -> (start, end)
    mode -> 
        'pass' -> Allow through only frequencies in band
        'block' -> Block filters in band
    
    plotband => (start,end) in frequencies for plots
    
    """
    nf = f.size
    
    #Pad
    fin = np.pad(f, pad_width=nf, mode='wrap')
    #Put a hamming window over the whole padded array
    fin = np.hanning(3*nf)*fin
    
    fft = np.fft.fft(fin)
    
    freq = np.fft.fftfreq(3*nf, d=dt)
    pfreq = freq[0:int(1.5*nf)]
    
    if band[0] is None:
        a = 0
    else:
        a = np.argmin(np.abs(band[0] - pfreq))
        
    if band[1] is None:
        b = int(1.5*nf)
    else:
        b = np.argmin(np.abs(band[1] - pfreq))
       

    mask = np.zeros(3*nf) 
    mask[a:b] =  np.hanning(b-a)
    mask[3*nf-b:3*nf-a] =  np.hanning(b-a)
    
    
    if mode == 'block':
        mask = np.ones(mask.size) - mask

         
        
    if plots:
        fig, ax = plt.subplots()
        ax.plot(freq, mask/np.max(mask))
        ax.plot(freq, np.abs(fft)/np.max(np.abs(fft)))
        
        ax.axvline(x=freq[a], color='green')
        ax.axvline(x=freq[b], color='red')
        
        if plotband is not None:
            ax.set_xlim(plotband )
        
        
        plt.show()
        
    fout = np.real(np.fft.ifft(fft*mask))
    #Divide back out the mask that was originally applied in time space
    
    fout = fout/np.hanning(3*nf)
    fout = fout[nf:2*nf] #Trim
    
    if plots:
        plt.plot(f)
        plt.plot(fout)
        plt.show()
    
    return fout
    
    

def lowpassFilter2D(arr, dx, dy, cutoff=None, order=6, plots=False):

    nx, ny = arr.shape
    fft = np.fft.fft2(arr)
    fft = np.fft.fftshift(fft)

    xfreq = np.fft.fftshift(np.fft.fftfreq(nx, d=dx))
    yfreq = np.fft.fftshift(np.fft.fftfreq(ny, d=dy))
    

    xfreq, yfreq = np.meshgrid(xfreq, yfreq, indexing='ij')
    
    
    freq = np.sqrt(np.power(xfreq,2) + np.power(yfreq,2))

    
    #Create a 2D butterworth filter
    filt = 1/np.sqrt(1 + np.power(freq/cutoff,2*order))
    
    if plots:
         plt.pcolormesh(filt)
         plt.show()
    
    #Apply to the data, then reverse each direction and apply again to
    #make the filter zero-phase
    fft = fft*filt
    fft = np.flip(fft, axis=0)
    fft = np.flip(fft, axis=1)
    fft = fft*filt
    fft = np.flip(fft, axis=0)
    fft = np.flip(fft, axis=1)
    

    fft = np.fft.ifftshift(fft)
    arr = np.real(np.fft.ifft2(fft))
    
    return arr  
    

    
    
    
    
if __name__ == '__main__':
    
    
    dk, x, y, arr = synthdata.wavey2D()
    
    #lowpassFilter2D(arr.T, dk, dk, cutoff=10, order=5, plots=True)
    
    f = arr[:,0]
    fftFilter(f, .1, plots=True, band=(1, 10))
    
    """
    dk, x, y, arr = synthdata.wavey2D()
    
    fig, ax = plt.subplots(figsize = [4, 4])
    cplot = ax.contourf(x, y, arr.T, levels=50, vmin=-1, vmax=1)
    
    
    arr = fftFilter(arr, dk, band=(15, None), mode='pass', axis=0, plots=False)
    
    fig, ax = plt.subplots(figsize = [4, 4])
    cplot = ax.contourf(x, y, arr.T, levels=50, vmin=-1, vmax=1)
    
    
    
    dt, t, arr = synthdata.twoWavePackets()
    
    f = fftFilter(arr, dt, band=(4e6, 8e6), mode='none', axis=0, plots=True)
    
    plt.plot(t*1e6, f)
    plt.xlim(0,100)
    plt.show()
    """

