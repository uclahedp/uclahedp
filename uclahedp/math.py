#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 10:57:18 2019

@author: peter
"""
import numpy as np

import h5py

def curl(arr, xax, yax, zax, xaxis, yaxis, zaxis):
    nx = xaxis.shape[0]
    ny = yaxis.shape[0]
    nz = zaxis.shape[0]


    #create output array, initially fill with NaN
    c = np.empty(  arr.shape  )
    c.fill(np.nan)
    
    if ny > 2 and nz > 2:
        print('Compute curl x component')
        dy = np.mean(np.gradient(yaxis))
        dz = np.mean(np.gradient(zaxis))
        c[..., 0] = np.gradient(c[..., 2], dy, axis=yax) - np.gradient(c[..., 1], dz, axis=zax)
    
    if nx > 2 and nz > 2:
        print('Compute curl y component')
        dx = np.mean(np.gradient(xaxis))
        dz = np.mean(np.gradient(zaxis))
        c[..., 1] = np.gradient(c[..., 0], dz, axis=zax) - np.gradient(c[..., 2], dx, axis=xax)

        
    if nx > 2 and ny > 2:
        print('Compute curl z component')
        dx = np.mean(np.gradient(xaxis))
        dy = np.mean(np.gradient(yaxis))
        c[..., 2] = np.gradient(c[..., 1], dx, axis=xax) - np.gradient(c[..., 0], dy, axis=yax)
        
    return c
        
        
        

if __name__ == '__main__':
    f = r'/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run103_PL11B_full.hdf5'
    
    with h5py.File(f, 'r') as sf:
        arr = sf['data'][0:10, ...]
        dimlabels = sf['data'].attrs['dimensions']
        xaxis = sf['xaxes'][:]
        yaxis = sf['yaxes'][:]
        zaxis = sf['zaxes'][:]

        xax = 1
        yax = 2
        zax= 3
        
        c = curl(arr, xax, yax, zax, xaxis, yaxis, zaxis)