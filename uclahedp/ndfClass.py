#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:05:25 2018

@author: peter
"""
import numpy as np
import h5py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:08:18 2018

@author: peter

parent: ndf 


"""

class ndf:
    """Parent class for a generic HEDP dataset of any dimensionality"""
    def __init__(self, filename): # Initialize object
        self.filename = filename
        self.f = h5py.File(filename, 'r')
        
        # Init quantities that every HDF file should have
        self.drive = self.f.attrs['drive']

    
class ndfarr (ndf):
    """
    Class of gridded NDF datasets
    Inherits ndf
    """
    def __init__(self, filename):
        ndf.__init__(self, filename)
        self.dim_names = self.f.attrs['pos_labels']
    
    
    
        
        
if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    obj = ndfarr(fname)
    print(type(obj))
    print(obj.dim_names)
    print(obj.drive)


        #    self.attrs = {}
         #   keys = list(f.attrs)
          #  for key in keys:
           #     self.attrs[key] = f.attrs[key]