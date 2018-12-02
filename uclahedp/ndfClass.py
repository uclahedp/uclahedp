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
from astropy import units as u
from collections import OrderedDict

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
        self.data = None
        self.axes = {}
        
        self.dimlabels = []
        self.dimunits = []
        self.data_label = None
        self.data_unit = None
        


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
        
        
        
        data_unit = f.attrs['data_unit']
        self.data =  f['data'][:] *  u.Unit(data_unit, parse_strict = 'warn' )
        self.data_label = f.attrs['data_label']


        dimlabels =  [ x.decode("utf-8") for x in f.attrs['dimlabels']  ]
        dimunits =  [ x.decode("utf-8") for x in f.attrs['dimunits']  ]
        self.ndim = len(dimlabels)
        
        self.axes = OrderedDict()
        for di in range(self.ndim):
            self.axes[ dimlabels[di] ] = ( f[ 'ax' + str(di) ][:] * u.Unit(dimunits[di], parse_strict = 'warn' ) )

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
            
        dimlabels = self.axes.keys()

        f.attrs['dimlabels'] = [s.encode('utf-8') for s in dimlabels] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        f.attrs['data_label'] = self.data_label
        f.attrs['data_unit'] = str(self.data.unit)

        
        dimunits = []
        for i, key in enumerate(dimlabels):
            f['ax' + str(i)] = self.axes[key].value
            dimunits.append( str( self.axes[key].unit ) )

        f.attrs['dimunits'] = [s.encode('utf-8') for s in dimunits] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Saving NDF object as HDF ' + self.save_filepath
        self.log.append(entry)
        f['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    
    
    def getAxis(self, key):
        key= self._getAxisKey(key)
        return self.axes[key]
    
    
    def avgDim(self, dimkey ):
        key = self._getAxisKey(dimkey)
        ax_ind = self._getAxisInd(dimkey)
   
        #Average the data
        self.data = np.average(self.data, axis=ax_ind)
        #Remove the axis from the the axes dictionary
        self.axes.pop(key)
        self.ndim = self.ndim-1
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Averaged over ' + key + ' axis'
        self.log.append(entry)
        



    def collapseDim(self, dimkey, value):
        #Find the axis index cooresponding to dim, and the axes it cooresponds to
        dimkey = self._getAxisKey(dimkey)
        ax_ind = self._getAxisInd(dimkey)
        ax = self.axes[dimkey]
        
        #If value isn't already an astropy quantity, assume it has the same units as the axis.
        if not isinstance(value, u.Quantity ):
            value = value * ax.unit

        #Find the index along the axis that is closes to 'value'
        ind = np.where(abs(ax - value) == abs(ax - value).min())[0][0]

        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Collapsed ' + dimkey + ' axis to ' + dimkey +  ' = ' + str(ax[ind])
        self.log.append(entry)


        #Collapse the dimension in the array
        self.data = np.take(self.data, ind, axis=ax_ind)
        #Delete the axis from the object
        self.axes.pop(dimkey)
        self.ndim = self.ndim-1
        
        
        
    def convertAxisUnit(self, dimkey, unit):
        #Coerce the key to be valid
        dimkey = self._getAxisKey(dimkey)
        
        
        
        if isinstance(unit, str ):
            try:
                unit = u.Unit( unit )
            except ValueError:
                print("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
                return None

        if not isinstance(unit, u.UnitBase):
            print("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
            return None
        
        
        old_unit = self.axes[dimkey].unit    
        self.axes[dimkey] = self.axes[dimkey].to(unit)


        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Converted ' + dimkey + ' from ' + str(old_unit) + ' to ' + str(self.axes[dimkey].unit) 
        self.log.append(entry)
        
        
        




    def plot(self, xrange=[None,None], yrange=[None,None], zrange=[None,None]):
        keyiter = iter( self.axes.keys() )
        if self.ndim == 1:
            print("Call simple 1D plotting routine")
            xkey = next(keyiter)
            
            #Convert axis range units
            if isinstance(xrange[0], u.Quantity):
                xrange[0] = xrange[0].to_value(self.axes[xkey].unit)
            if isinstance(xrange[1], u.Quantity):
                xrange[1] = xrange[1].to_value(self.axes[xkey].unit)
            if isinstance(yrange[0], u.Quantity):
                yrange[0] = yrange[0].to_value(self.data.unit)
            if isinstance(yrange[1], u.Quantity):
                yrange[1] = yrange[1].to_value(self.data.unit)
       

            plt.plot(self.getAxis(xkey), self.data)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(xkey + ' (' + str(self.axes[xkey].unit) + ')'  )
            plt.ylabel(self.data_label.title() + ' (' + str(self.data.unit) + ')'  )
            plt.show()
        elif self.ndim == 2:
            print("Call simple 2D contour plotting routine")
            xkey = next(keyiter)
            ykey = next(keyiter)
            
        elif self.ndim == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
    
    
    #**************
    #HELPER METHODS
    #**************
    def _getAxisKey(self, dimkey):
        key = [k for k in  self.axes.keys() if k.lower() == dimkey.lower()]
        if len(key) == 0:
            print("No axis exists with the name: " + dimkey)
            return None
        elif len(key) > 1:
            print("Multiple axes exist with the name:" + dimkey)
            return None
        else:
            return key[0]
        
        
    def _getAxisInd(self, key):
        ind = [i for i,k in  enumerate(self.axes.keys()) if k.lower() == key.lower()]
        if len(ind) == 0:
            print("No axis exists with the name: " + key)
            return None
        elif len(ind) > 1:
            print("Multiple axes exist with the name:" + key)
            return None
        else:
            return ind[0]





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
    
    print(np.shape(obj.data))
    obj.avgDim('Reps')
    print(np.shape(obj.data))
    
    
    #print(dir(obj))
    print(np.shape(obj.getAxis('X') ))
    
    print(np.shape(obj.data))
    obj.collapseDim('x', .05*u.m)
    print(np.shape(obj.data))
    obj.collapseDim('channels', 1)
    print(np.shape(obj.data))
    
    obj.convertAxisUnit('time', u.us)
 
    
    obj.plot(xrange=[0, 40], yrange=[-150*u.mV, 150*u.mV])
    
    obj.saveHDF(sname)


    