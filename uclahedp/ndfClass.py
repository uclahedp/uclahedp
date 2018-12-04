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
        self.data_label = None
        self.log = []

        


    def readHDF(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:   
            self.unpack(f)
        
    def saveHDF(self, filepath):
        self.save_filepath = filepath
        with h5py.File(filepath, 'w') as f:
            self.pack(f)


        
    def unpack(self, f):
        self._log('Opened HDF file as NDF object: ' + self.read_filepath)

        self.data =  f['data'][:] *  u.Unit(f.attrs['data_unit'], parse_strict = 'warn' )
        self.data_label = f.attrs['data_label']


        dimlabels =  [ x.decode("utf-8") for x in f.attrs['dimlabels']  ]
        dimunits =  [ x.decode("utf-8") for x in f.attrs['dimunits']  ]
        self.axes = []
        for i in range(len(dimlabels)):
            ax = f[ 'ax' + str(i) ][:] * u.Unit(dimunits[i], parse_strict = 'warn' )
            name = dimlabels[i]
            a = {'name':name, 'axis':ax}
            self.axes.append(a)

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


        f.attrs['dimlabels'] = [ax['name'].encode('utf-8') for ax in self.axes] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        f.attrs['data_label'] = self.data_label
        f.attrs['data_unit'] = str(self.data.unit)

        
        dimunits = []
        for i, ax in enumerate(self.axes):
            f['ax' + str(i)] = ax['axis'].value
            dimunits.append( str( ax['axis'].unit ) )

        f.attrs['dimunits'] = [s.encode('utf-8') for s in dimunits] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        self._log('Saving NDF object as HDF ' + self.save_filepath)
        f['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    
    
    def getAxis(self, name):
        name = str(name)
        return [ax['axis'] for ax in self.axes if ax['name'].lower().strip() == name.lower().strip()  ][0]
    
    def getAxisInd(self, name):
        name = str(name)
        return [i for i,ax in enumerate(self.axes) if ax['name'].lower().strip() == name.lower().strip()  ] [0]
    
    
    def avgDim(self, name ):
        ax_ind = self.getAxisInd(name)
        #Average the data
        self.data = np.average(self.data, axis=ax_ind)
        #Remove the axis from the the axes dictionary
        self.axes.pop(ax_ind)
        
        self._log('Averaged over ' + name + ' axis')
        
        



    def collapseDim(self, name, value):
        # Find the axis index cooresponding to dim, and the axes it cooresponds to
        ax_ind = self.getAxisInd(name)
        ax = self.getAxis(name)

        #If value isn't already an astropy quantity, assume it has the same units as the axis.
        if not isinstance(value, u.Quantity ):
            value = value * ax.unit

        #Find the index along the axis that is closes to 'value'
        ind = np.where(abs(ax - value) == abs(ax - value).min())[0][0]

        #Collapse the dimension in the array
        self.data = np.take(self.data, ind, axis=ax_ind)
        #Delete the axis from the object
        self.axes.pop(ax_ind)

        
        self._log('Collapsed ' + name + ' axis to ' + name +  ' = ' + str(ax[ind]))
        
        
        
        
        
    def thinDim(self, name, bin=10):
        # Get the index that goes with this name
        ax_ind = self.getAxisInd(name)
        ax = self.getAxis(name)
        # Get the current shape vector
        shape = list( self.data.shape )
        # Calculate the new length of dim for the given bin
        nbins  = int( shape[ax_ind] / bin )
        shape[ax_ind] = nbins
        # Construct a new shape vector and new data array
        arr = np.zeros(shape)
        newax = np.zeros(nbins)

        # Fill the new data array while doing the averaging
        for i in range(nbins):
            # Create a range of indices within this row to average
            indrange = range(i*bin, (i+1)*bin-1)
            # Take them out
            s = np.take(self.data, indrange, axis=ax_ind)
            # Average them
            s = np.average(s, axis=ax_ind)
            # Create a slice object for this slice
            pslice = [ slice(None) ]*arr.ndim
            pslice[ax_ind] = slice(i, (i+1), None)
            # Put the values into the slice
            arr[ tuple(pslice) ] = s
            
            #Now do the same thing for the cooresponding axis
            s = np.take(ax, indrange)
            s = np.average(s)
            np.put(newax, i, s)
        
 
            
        # Put the data back into the array
        self.data = arr*self.data.unit
        # Put the new axis back too
        self.axes[ax_ind]['axis'] = newax * self.axes[ax_ind]['axis'].unit
        
        
        
        
    def convertAxisUnit(self, name, unit):
        if isinstance(unit, str ):
            try:
                unit = u.Unit( unit )
            except ValueError:
                print("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
                return None

        if not isinstance(unit, u.UnitBase):
            print("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
            return None
        
        ax_ind = self.getAxisInd(name)
        ax = self.getAxis(name)
        
        old_unit = ax.unit    
        self.axes[ax_ind]['axis'] = ax.to(unit)

        self._log('Converted ' + name + ' from ' + str(old_unit) + ' to ' + str( (self.axes[ax_ind]['axis']).unit ) )
        
        
        




    def plot(self, xrange=[None,None], yrange=[None,None], zrange=[None,None]):
        xname = self.axes[0]['name']
        xaxes = self.axes[0]['axis']
        
        if len(self.axes) == 1:
            print("Call simple 1D plotting routine")
            
            
            #Convert axis range units
            if isinstance(xrange[0], u.Quantity):
                xrange[0] = xrange[0].to_value(self.axes[xkey].unit)
            if isinstance(xrange[1], u.Quantity):
                xrange[1] = xrange[1].to_value(self.axes[xkey].unit)
            if isinstance(yrange[0], u.Quantity):
                yrange[0] = yrange[0].to_value(self.data.unit)
            if isinstance(yrange[1], u.Quantity):
                yrange[1] = yrange[1].to_value(self.data.unit)
       
            #Make plot
            plt.plot(xaxes, self.data)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(xname + ' (' + str(xaxes.unit) + ')'  )
            plt.ylabel(self.data_label.title() + ' (' + str(self.data.unit) + ')'  )
            plt.show()
        elif len(self.axes) == 2:
            print("Call simple 2D contour plotting routine")
            yname = self.axes[1]['name']
            yaxes = self.axes[1]['axis']
            
            #Convert axis range units
            if isinstance(xrange[0], u.Quantity):
                xrange[0] = xrange[0].to_value(self.axes[xkey].unit)
            if isinstance(xrange[1], u.Quantity):
                xrange[1] = xrange[1].to_value(self.axes[xkey].unit)
            if isinstance(yrange[0], u.Quantity):
                yrange[0] = yrange[0].to_value(self.axes[ykey].unit)
            if isinstance(yrange[1], u.Quantity):
                yrange[1] = yrange[1].to_value(self.axes[ykey].unit)
            if isinstance(zrange[0], u.Quantity):
                zrange[0] = zrange[0].to_value(self.data.unit)
            if isinstance(zrange[1], u.Quantity):
                zrange[1] = zrange[1].to_value(self.data.unit)
            
            print(xkey)
            print(np.shape(self.getAxis(xkey)))
            print(ykey)
            print(np.shape(self.getAxis(ykey)))
            print(np.shape(self.data))
            plt.figure()
            plt.contourf(xaxes, yaxes, self.data.T)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(xname + ' (' + str(xaxes.unit) + ')'  )
            plt.ylabel(yname + ' (' + str(yaxes.unit) + ')'  )
            plt.title(self.data_label.title() + ' (' + str(self.data.unit) + ')'  )
            plt.show()
            
        elif len(self.axes) == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
    
    
    #**************
    #HELPER METHODS
    #**************
    def _log(self, message):
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + str(message)
        self.log.append(entry)
        





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
    sname = r"C:\Users\Peter\Desktop\TempData\testsave.h5"

    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    fname = r"C:\Users\Peter\Desktop\TempData\run56_LAPD1_pos_raw.h5"
    #fname = r"C:\Users\Peter\Desktop\TempData\run102_PL11B_pos_raw.h5"
    
    obj = ndf()
    
    obj.readHDF(fname)
    #obj.plot()
    #obj.data = obj.data*0 + 20 #Put in some fake data to test that it changes
    
    #print(obj.getAxis('reps'))
    
    #obj.thinDim('time', bin=20)
    
    print(np.shape(obj.data))
    obj.avgDim('Reps')
    print(np.shape(obj.data))
    
    
    #print(dir(obj))
    #print( obj.getAxis('time') )
    
    print(np.shape(obj.data))
    obj.collapseDim('x', 5*u.cm)
    
    print(np.shape(obj.data))
    obj.collapseDim('channels', 2)
    print(np.shape(obj.data))
    
    obj.convertAxisUnit('time', u.us)

    
    print(np.shape(obj.data))
    
    
    obj.plot(xrange=[0,20 ])
    #obj.plot( xrange = [0, 40*u.us])
    
    obj.saveHDF(sname)


    