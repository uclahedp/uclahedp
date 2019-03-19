#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 11:45:28 2019

@author: peter
"""
import numpy as np
from scipy.signal import gaussian as gaussian
import matplotlib.pyplot as plt


def twoWavePackets(verbose=False, plots=False):
    # Create a test signal with two frequency components
    n = int(1e6)
    f1, f2 = 6e6,2e6
    dt=1e-10
    t = np.arange(n)*dt
    loc1 = np.argmin(np.abs(t*1e6-20))
    width1 = int(20*1e-6/dt/2.)
    loc2 = np.argmin(np.abs(t*1e6-40))
    width2 = int(80*1e-6/dt/2.)
    m1 = np.zeros(n)
    m1[loc1-width1:loc1+width1] = gaussian(width1*2, width1/12)
    m2 = np.zeros(n)
    m2[loc2-width2:loc2+width2] = gaussian(width2*2, width2/12)

    arr = np.cos(2*np.pi*f1*t)*m1  + 2*np.cos(2*np.pi*f2*t)*m2
    
    if plots:
        plt.plot(t* 1e6, arr)
        plt.show()
    
    return dt, t, arr


def cpWave(plots=False, rcp=None):
    
    if rcp is None or rcp is True:
        rcp = True
    else:
        rcp = False #Generate LCP
    
    n = int(1e5)
    f1 = 5e5
    dt=1e-10
    t = np.arange(n)*dt
    
    if rcp:
        ax = np.sin(2*np.pi*f1*t)
        ay = np.cos(2*np.pi*f1*t)
    else:
        ay = np.sin(2*np.pi*f1*t)
        ax = np.cos(2*np.pi*f1*t)
        
    if plots:
        plt.plot(t*1e6, ax, '-', t*1e6, ay, '--')
        
    return t, ax, ay


if __name__=='__main__':
    #x = twoWavePackets(plots=True)
    x = rcpWave(plots=True)