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
import h5py

import os

def polarization_decomp(bx,by, flipz=False):
    """
    INPUTS
    bx, by (Float arrays)
    Transverse field components
    
    flipz (Boolean)
    If False, wave is assumed to be propagating along the z direction defined
    by the given x and y directions and the right-hand rule. If True, the wave
    is assumed to propagate in the -Z direction. 
    
    OUPTUTS
    br, bl (Float arrays)
    Polarization decomposition of the field: decomposed onto RCP and LCP
    basis vectors.
    """

    # Compute forward Fourier transforms
    bxf = fft.fft(bx)
    byf = fft.fft(by)
    # Get phase-shifted superpositions (byf+I*bzf)/2,(byf-I*bzf)/2
    if flipz:
        brf = (0.5+0.j)*bxf + (0.+0.5j)*byf
        blf = (0.5+0.j)*bxf - (0.+0.5j)*byf
    else:
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
    t, ax, ay = cpWave(rcp=True)
    br, bl = polarization_decomp(ax, ay)
    
    fig, ax = plt.subplots()
    ax.plot(t*1e6, br, '-', t*1e6, bl, '--')
    

def hodograph(ax, ay, indices=None):
    if indices is None:
        indices = np.arange(ax.shape[0]-1)
        
        
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
    t, ax, ay = cpWave(rcp=True)
    
    indices = np.arange(0, 15000, 500)
    x,y,u,v = hodograph(ax, ay, indices=indices)
    
    fig, ax = plt.subplots()
    ax.set_aspect(1.0)
    
    for i in range(x.size):
        ax.arrow(x[i], y[i], u[i], v[i])
    
    #ax.quiver(x,y,u,v, pivot='tail')





if __name__ == '__main__':
    #_test_polarization_decomp()
   _test_hodograph()
   
   
   #Test Hodograph
   """
   f = os.path.join("G:", "LAPD_Mar2018", "FULL", "run40_LAPD7_full.hdf5")
   
   times = np.arange(5, 40, 0.5)
   tinds = np.zeros([times.size], dtype=int)
   
   with h5py.File(f) as f:
        bt = f['time'][:]*1e6
        bx = f['data'][0,:,0]
        by = f['data'][0,:,1]
        
   for i,t in enumerate(times):
       tinds[i] = np.argmin(np.abs(bt - t))
       
        
   x0,y0,dx,dy = hodograph(bx, by, indices=tinds)
        
   print(x0.shape)
   
   fig, ax = plt.subplots()
   
   for i in range(times.size):
       ax.arrow(x0[i], y0[i], dx[i], dy[i])
   """
    
   """
    #THIS IS ALL FOR TESTING DECOMP
    title='Parallel'
    f = os.path.join("G:", "LAPD_Mar2018", "FULL", "run40_LAPD7_full.hdf5")
    #f =  '/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run56_LAPD1_full.hdf5'
    
    
    #title='Perpendicular'
    #f = '/Volumes/PVH_DATA/LAPD_Jan2019/FULL/run34_LAPD10_full.hdf5'
    
    with h5py.File(f) as f:
        
        t = f['time'][:]
        
        bx = f['data'][0,:,0]
        by = f['data'][0,:,1]
        #bx = f['data'][:,30,0,0,0,0]
        #by = f['data'][:,30,0,0,0,1]
        

        fig, ax = plt.subplots()
        br, bl = polarization_decomp(bx, by, flipz=True)
        ax.plot(t*1e6, br, 'b-', label='RCP' )
        ax.plot(t*1e6, bl, 'r-', label='LCP')
        ax.set_xlim((0,50))
        ax.set_xlabel('t (us)')
        ax.set_ylabel('dB (G)')
        ax.set_title(title)
        ax.legend(loc=4)

    """ 
        
        
   
        