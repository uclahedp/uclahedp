#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:05:25 2018

@author: peter
"""
import numpy as np
import h5py
import datetime
import matplotlib.pyplot as plt   

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:08:18 2018

@author: peter

parent: ndf 


"""

class ndf:
    """Parent class for a generic HEDP dataset of any form"""

    def __init__(self): # Initialize object
        #TODO: Initialize all of the variables here for use in a separate constructor method?
        pass


    def readHDF(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:   
            self.unpack(f)
        
    def saveHDF(self, filepath):
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
        ind = self._getAxisInd(dim)
        return self.axes[ind]
    
    
    def avgDim(self, dim ):
        ind = self._getAxisInd(dim)
        print(np.shape(self.data))
        print(ind)
        #Average the data
        self.data = np.average(self.data, axis=ind)
        
        self._deleteAxisRefs(ind)
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Averaged over dimension ' + str(ind) + ' (' + self.dimlabels[ind] + ')' 
        self.log.append(entry)
        



    def collapseDim(self, dim, value):
        #Find the axis index cooresponding to dim, and the axes it cooresponds to
        ax_ind = self._getAxisInd(dim)
        ax = self.axes[ax_ind]
        
        #Find the index along the axis that is closes to 'value'
        ind = np.where(abs(ax - value) == abs(ax - value).min())[0][0]
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Collapsed dimension ' + str(ax_ind) + ' (' + self.dimlabels[ax_ind] + ')' 
        self.log.append(entry)

        #Collapse the dimension in the array
        self.data = np.take(self.data, ind, axis=ax_ind)
        #Delete the axis from the object
        self._deleteAxisRefs(ax_ind)
        




    def plot(self, xrange=(None,None), yrange=(None,None), zrange=(None,None)):
        if self.ndim == 1:
            print("Call simple 1D plotting routine")
            plt.plot(self.getAxis(0), self.data)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(self.dimlabels[0] + ' (' + self.units[0] + ')'  )
            #plt.ylabel('BY (G)')
            plt.show()
        elif self.ndim == 2:
            print("Call simple 2D contour plotting routine")
            
        elif self.ndim == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
    
    
    #**************
    #HELPER METHODS
    #**************
    def _getAxisInd(self, dim):
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



    def _deleteAxisRefs(self, ind):
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
    
    
    obj = ndf()
    obj.readHDF(fname)
    #obj.plot()
    #obj.data = obj.data*0 + 20 #Put in some fake data to test that it changes
    
    obj.avgDim('Reps')
    #print(dir(obj))
    #print(np.shape(obj.getAxis('X') ))
    obj.collapseDim('x', 1)
    obj.collapseDim('channels', 1)
    
    obj.plot(xrange=(0,40*1e-6))
    
    obj.saveHDF(sname)


    