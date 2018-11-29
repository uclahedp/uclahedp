#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:05:25 2018

@author: peter
"""
import numpy as np
import h5py

class ndf:
    
    def __init__(self, f):
        with h5py.File(rawfilepath, 'r') as f:
            
            self.gridded = f.attrs['gridded']
            self.data_form = f.attrs['data_form']
            
            self.data = f['data'][:]
            
            self.t0 = f['t0'][:]
            
            if self.gridded:
                self.xgv = f['grid_xgv'] # Change to just xgv in hdf?
                self.ygv = f['grid_ygv']
                self.zgv = f['grid_zgv']
                
            if self.data_form == 'pos':
                self.pos = f['pos'][:]
            
            
            self.data_unit = f['data'].attrs['unit']
            self.pos_unit = f['pos'].attrs['unit']
            self.t0_unit = f['t0'].attrs['unit']
            
            
            self.attrs = {}
            keys = list(f.attrs)
            for key in keys:
                self.attrs[key] = f.attrs[key]
                
    
    def save(self, filepath):
        with h5py.File(filepath, 'w') as f:
            print(1)
            

            
            

        
        
if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    
    obj = ndf(fname)
    
    print(type(obj.data))
    print(np.shape(obj.data))