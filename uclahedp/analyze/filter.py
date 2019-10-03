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



def fftFilter(f, dt, band=(None,None), axis=0, mode='pass', plots=False):
    """
    band -> (start, end)
    mode -> 
        'pass' -> Allow through only frequencies in band
        'block' -> Block filters in band
    """

    
    nf = f.shape[axis]
    fft = np.fft.fft(f, axis=axis)
    freq = np.fft.fftfreq(nf, d=dt)
    
    
    #Find endpoints of the bandpass window
    if band[0] is None or band[0] ==0 :
        a = 1
    else:
        a = np.argmin(np.abs(freq - band[0]))
           
    if band[1] is None:
        b = int(nf/2) #Set to the nyquist frequency (highest freq bin)
    else:
        b = np.argmin(np.abs(freq - band[1]))

    mask = np.zeros(fft.shape)
    
    
    #Assemble mask slice
    posslice = []
    negslice = []
    
    maskshape = []
    for i,n in enumerate(mask.shape):
        if i==axis:
            posslice.append( slice(a,b,None) )
            negslice.append( slice(-b,-a,None) )
            maskshape.append(b-a)
        else:
            posslice.append( slice(None,None,None) )
            negslice.append( slice(None,None,None) )
            maskshape.append(1)
    mask[tuple(posslice)] = np.hanning(b-a).reshape(maskshape)
    mask[tuple(negslice)] = np.hanning(b-a).reshape(maskshape)
    
    
    
    if plots:
        fig, ax = plt.subplots(figsize = [4, 4])
        if fft.ndim == 1:
            ax.plot(freq, fft/np.max(fft))
            ax.plot(freq, mask)
            ax.set_xlim(0.25*freq[a],4*freq[b])
    
    if mode == 'pass':
        fft = fft*mask
    elif mode == 'block':
        fft = fft - fft*mask
    elif mode == 'none':
        #Including this case for testing
        pass
    

        
    f = np.fft.ifft(fft, axis=axis)

    return f
        
   
   
    

def lowpassFilter2D(arr, dx, dy, cutoff=10, plots=False):

    
    nx, ny = arr.shape
    fft = np.fft.fft2(arr)
    fft = np.fft.fftshift(fft)
    
    
    xfreq = np.fft.fftfreq(nx, d=dx)[0:int(nx/2.0)]
    yfreq = np.fft.fftfreq(ny, d=dy)[0:int(ny/2.0)]
    
    xi = np.argmin(np.abs(xfreq - cutoff))
    yi = np.argmin(np.abs(yfreq - cutoff))
    
    a = xfreq[0:xi].size*2
    b = yfreq[0:yi].size*2
    
    #Deal with odd/even sizes here so padding will be integers
    if nx % 2 == 1:
        a = a - 1 
    if ny % 2 == 1:
        b = b - 1 
    
    mask = np.outer(np.hanning(a), np.hanning(b))
    
    if plots:
        plt.pcolormesh(mask)
        plt.show()
    
    xpad = int((nx - a)/ 2)
    ypad = int((ny- b)/2)
    mask = np.pad(mask, pad_width=((xpad,xpad), (ypad,ypad)), mode='constant')
    
    fft = np.fft.ifftshift(mask*fft)
    arr = np.fft.ifft2(fft)


    return arr


    
    

    
    
    
    
if __name__ == '__main__':
    
    
    dk, x, y, arr = synthdata.wavey2D()
    
    lowpassFilter2D(arr, dk, dk, cutoff = 15)
    
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

