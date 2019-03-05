#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 18:02:45 2019

@author: Martin Weidl and Peter Heuer
"""

# Python module imports
import numpy as np
import scipy.fftpack as fft
from uclahedp.tests.synthdata import cpWave as cpWave
import matplotlib.pyplot as plt

def polarization_decomp(bx,by):
    """
    IN:  by, bz Real 1D-arrays with perpendicular field components
    OUT: br, bl Complex 1D-arrays with the right- and left-hand circularly 
                polarized components of by, bz, with respect to the 
                positive x-axis. The real parts correspond to by and the 
                imaginary parts to bz, which should only differ by much if 
                the signal is strongly elliptically polarized.
    """
    # Compute forward Fourier transforms
    bxf = fft.fft(bx)
    byf = fft.fft(by)
    # Get phase-shifted superpositions (byf+I*bzf)/2,(byf-I*bzf)/2
    blf = (0.5+0.j)*bxf + (0.+0.5j)*byf
    brf = (0.5+0.j)*bxf - (0.+0.5j)*byf
    # Store length of arrays in array_len (assumed to be even)
    array_len = int( np.size(bxf) /2. )
    
    print(array_len/2)
    # Set all negative frequencies (equiv to n>=N/2) to zero
    brf[array_len:] = 0
    blf[array_len:] = 0
    # Set the zero-frequency component to zero
    brf[0] = 0
    blf[0] = 0
    # Compute inverse Fourier transforms
    br = fft.ifft(brf)
    bl = fft.ifft(blf)
    # Return (br,bl)
    return (br,bl)


def _test_polarization_decomp():
    t, ax, ay = cpWave(rcp=False)
    br, bl = polarization_decomp(ax, ay)
    plt.plot(t*1e6, br, '-', t*1e6, bl, '--')
    

def hodograph(ax, ay, indices=None):
    if indices is None:
        indices = np.arange(ax.shape[0])
        
        
    x0 = np.zeros(indices.shape[0])
    y0 = np.zeros(indices.shape[0])
    dx = np.zeros(indices.shape[0])
    dy = np.zeros(indices.shape[0])

    
    for i, ind in enumerate(indices):

        x0[i] = ax[ind]
        y0[i] = ay[ind]
            
        dx[i] = ax[ind+1] - x0[i]
        dy[i] = ay[ind+1] - y0[i]
        
    return x0,y0,dx,dy
    
def _test_hodograph():
    t, ax, ay = cpWave(rcp=False)
    indices = np.arange(0, 15000, 500)
    x,y,u,v = hodograph(ax, ay, indices=indices)
    
    plt.quiver(x,y,u,v, pivot='tail')
    plt.show()




if __name__ == '__main__':
    _test_polarization_decomp()
    _test_hodograph()