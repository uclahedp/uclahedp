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
    def __init__(self, f): # Initialize object

        # Init quantities that every HDF file should have
        self.drive = f.attrs['drive']
        # TODO others here...
        

    
class ndfgrid (ndf):
    """
    Class of gridded NDF datasets
    Inherits ndf
    """
    def __init__(self, f):
        ndf.__init__(self, f)
        self.dimlabels =  [ x.decode("utf-8") for x in f.attrs['dimlabels']  ]
        self.units =  [ x.decode("utf-8") for x in f.attrs['units']  ]
        
        self.ndim = len(self.dimlabels)
        
        #Add quantities we only want if a dimension is non-trivial
        for dim in self.dimlabels:
            if dim is 'time':
                self.time = f['time']
                self.nti = len(self.time)
            if dim is 'x':
                self.x = f['x']
                self.nx = len(self.x)
            if dim is 'y':
                self.y = f['y']
                self.ny = len(self.y)
            if dim is 'z':
                self.x = f['z']
                self.nz = len(self.z)
            if dim is 'reps':
                self.reps = f['reps']
            if dim is 'channels':
                self.channels = f['channels']
        
        
        
class ndfpoints(ndf):
    """
    Class of NDF dataset that is not gridded
    Inherits ndf
    """
    def __init__(self, f):
        pass
    
        
        
if __name__ == "__main__":
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    f = h5py.File(fname, 'r')
    obj = ndfgrid(f)
    f.close()
    print(type(obj))
    print(obj.dimlabels)
    print(obj.drive)


        #    self.attrs = {}
         #   keys = list(f.attrs)
          #  for key in keys:
           #     self.attrs[key] = f.attrs[key]