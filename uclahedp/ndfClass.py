#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:05:25 2018

@author: peter
"""
import numpy as np
import h5py
import datetime

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:08:18 2018

@author: peter

parent: ndf 


"""

class ndf:
    """Parent class for a generic HEDP dataset of any dimensionality"""

    def __init__(self, filepath): # Initialize object
        self.read(filepath)


    def read(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:   
            self.unpack(f)
        
    def save(self, filepath):
        self.save_filepath = filepath
        with h5py.File(filepath, 'w') as f:
            elf.pack(f)


        
    def unpack(self, f):
        self.log =  [ x.decode("utf-8") for x in f['log'] ]
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Opened HDF file as NDF object: ' + self.read_filepath
        self.log.append(entry)
        

        self.data = f['data'][:]

        self.dimlabels =  [ x.decode("utf-8") for x in f.attrs['dimlabels']  ]
        self.units =  [ x.decode("utf-8") for x in f.attrs['units']  ]
        self.ndim = len(self.dimlabels)
        
        # Init quantites that we want to have easily available
        self.run = f.attrs['run']
        self.drive = f.attrs['drive']
        self.daq = f.attrs['daq']
        self.dt = f.attrs['dt']
        self.probe_name = f.attrs['probe_name']
        self.probe_type = f.attrs['probe_type']
        # TODO others here...
        

        #Make a list of all of the attributes in the HDF, so we can perpetuate them
        self.attrs = {}
        keys = list(f.attrs)
        for key in keys:
            self.attrs[key] = f.attrs[key]
            
        
    def pack(self, f):
        f['data'] = self.data

        # Add the original list of all attributes back to the dataset
        # This comes before the explictly coded ones, so those will take 
        # precedence in case of changes
        for key in self.attrs.keys():
            f.attrs[key] =self.attrs[key]
        
        # Write the important attributes, with the new object values overwriting
        # any duplicates from self.attrs
        f.attrs['run'] = self.run
        f.attrs['drive'] = self.drive
        f.attrs['daq'] =  self.daq 
        f.attrs['dt'] =  self.dt
        f.attrs['probe_name'] = self.probe_name 
        f.attrs['probe_type'] = self.probe_type
        
        f.attrs['dimlabels'] = [s.encode('utf-8') for s in self.dimlabels] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        f.attrs['units'] = [s.encode('utf-8') for s in self.units] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Saving NDF object as HDF ' + self.save_filepath
        self.log.append(entry)
        f['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    


    



class ndfgrid (ndf):
    """
    Class of gridded NDF datasets
    Inherits ndf
    """
    def __init__(self, filepath):
        ndf.__init__(self, filepath)
        
        
    def read(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:
            self.unpack(f)
        
    def save(self, filepath):
        self.save_filepath = filepath
        with h5py.File(filepath, 'w') as f:
            self.pack(f)
        
      
        
    def unpack(self,f):
        ndf.unpack(self,f)
        #Add quantities we only want if a dimension is non-trivial
        if 'time' in self.dimlabels:
            self.time = f['time'][:]
            self.nti = len(self.time)
        else:
            self.nti = 1
        if 'x' in self.dimlabels:
            self.x = f['x'][:]
            self.nx = len(self.x)
        else: 
            self.nx = 1
        if 'y' in self.dimlabels:
            self.y = f['y'][:]
            self.ny = len(self.y)
        else:
            self.ny = 1
        if 'z' in self.dimlabels:
            self.z = f['z'][:]
            self.nz = len(self.z)
        else:
            self.nz = 1
        if 'reps' in self.dimlabels:
            self.reps = f['reps'][:]
            self.nreps = len(self.reps)
        else:
            self.nreps = 1
        if 'channels' in self.dimlabels:
            self.channels = f['channels'][:]
            self.nchan = len(self.channels)
        else:
            self.nchan = 1
        self.npos = self.nx*self.ny*self.nz
        self.nshots = self.npos*self.nreps
               
        
        
        
        
    def pack(self, f):
        ndf.pack(self,f)
        #Add quantities we only want if a dimension is non-trivial
        for dim in self.dimlabels:
            print(dim)
            if dim == 'time':
                f['time'] = self.time
            if dim == 'x':
                f['x'] = self.x
            if dim == 'y':
                f['y'] = self.y 
            if dim == 'z':
                f['z'] = self.z 
            if dim == 'reps':
                f['reps'] = self.reps
            if dim == 'channels':
                f['channels'] = self.channels
           
            
            
            
    def plot(self):
        if self.ndim == 1:
            print("Call simple 1D plotting routine")
            
        elif self.ndim == 2:
            print("Call simple 2D contour plotting routine")
            
        elif self.ndim == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
        







        
class ndfpoints(ndf):
    """
    Class of NDF dataset that is not gridded
    Inherits ndf
    """
    def __init__(self, f):
        ndf.__init__(self, f)
        
        self.pos = f['pos'] # pos is only stored for non-gridded data which requires it.
        self.npos = len(self.pos)
        
        
        
    def save(self, f):
        ndf.save(self,f)
    
        
        
if __name__ == "__main__":
    sname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/testsave.h5"
    
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    
    obj = ndfgrid(fname)
    obj.plot()
    obj.save(sname)


    