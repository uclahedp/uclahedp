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
    """Parent class for a generic HEDP dataset of any form"""

    def __init__(self, filepath): # Initialize object
        self.read(filepath)


    def read(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:   
            self.unpack(f)
        
    def save(self, filepath):
        self.save_filepath = filepath
        with h5py.File(filepath, 'w') as f:
            self.pack(f)


        
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
        
        self.axes = []
        for di in range(len(self.dimlabels)):
            self.axes.append( f[ 'ax' + str(di) ][:] )

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
        
        
        for di in range(len(self.dimlabels)):
            f['ax' + str(di)] = self.axes[di]

        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Saving NDF object as HDF ' + self.save_filepath
        self.log.append(entry)
        f['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    
    
    def getAxis(self, dim):
        ind = self.getAxisInd(dim)
        return self.axes[ind]
    
    
    def avg(self, dim ):
        ind = self.getAxisInd(dim)
        #Average the data
        self.data = np.average(self.data, axis=ind)
        
        self.deleteAxisRefs(ind)
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Averaged over dimension ' + str(ind) + ' (' + self.dimlabels[ind] + ')' 
        self.log.append(entry)
        



    def collapseDim(self, dim, value):
        
        #TODO THIS CODE DOESN'T WORK
        
        ax_ind = self.getAxisInd(dim)
        ax = self.getAxis(dim)
        
        ax = self.axes[ax_ind]
        #ind = np.find_nearest(ax, value)
        ind = np.where(abs(ax - value) == abs(ax - value).min())[0][0]
        
        print(ind)
        print(ax[ind])

        print(np.shape(self.data))
        idx = tuple(  slice( s[0], s[1], None  ) for s in self.data )
        print(idx)
        self.data = self.data[s] #Must be a tuple for slicing
        #self.deleteAxisRefs(ax_ind)
        print(np.shape(self.data))
        
        
        
        
        
    
    def plot(self):
        if self.ndim == 1:
            print("Call simple 1D plotting routine")
            
        elif self.ndim == 2:
            print("Call simple 2D contour plotting routine")
            
        elif self.ndim == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
    
    
    #**************
    #HELPER METHODS
    #**************
    def getAxisInd(self, dim):
        if isinstance(dim, str):
            if not dim.isnumeric(): #Assume a dimlabel was requested
                ind = [i for i,v in enumerate(self.dimlabels) if v.lower() == dim.lower()]
                if len(ind) == 0:
                    print("No axis exists with the name: " + dim)
                    return None
                elif len(ind) > 1:
                    print("Multiple axes exist with the name:" + dim)
                    return None
                else:
                    ind=ind[0]
                        
        #If not a string, assume an integer index was given
        else:
            ind = int(dim)
            if ind >= self.ndim or ind < 0:
                print("Invalid dimension index to average: " + str(ind))
                return None
        return ind

    def deleteAxisRefs(self, ind):
        #Remove that axis from the relevant attribute lists
        self.axes.pop(ind)
        labelpopped = self.dimlabels.pop(ind)
        self.units.pop(ind)
        self.ndim = self.ndim - 1



    



class ndfgrid (ndf):
    """
    Class of gridded NDF datasets
    Inherits ndf
    """

    def unpack(self,f):
        ndf.unpack(self,f)
        
    def pack(self, f):
        ndf.pack(self,f)


        
class ndfpoints(ndf):
    """
    Class of NDF dataset that is not gridded
    Inherits ndf
    """
    def unpack(self, f):
        ndf.unpack(self,f)
        self.pos = f['pos'] # pos is only stored for non-gridded data which requires it.
        self.npos = len(self.pos)
        
        
    def pack(self, f):
        ndf.pack(self,f)
        
    
        
        
if __name__ == "__main__":
    sname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/testsave.h5"

    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    
    obj = ndf(fname)
    #obj.plot()
    obj.data = obj.data*0 + 20 #Put in some fake data to test that it changes
    
    #obj.avg('reps')
    #print(dir(obj))
    #print(np.shape(obj.getAxis('X') ))
    
    obj.collapseDim('X', 5)
    obj.save(sname)


    