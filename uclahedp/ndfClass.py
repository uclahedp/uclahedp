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
import copy
import hedpConstants as const

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:08:18 2018

@author: peter

parent: ndf 


"""

class ndf:
    """Parent class for a generic HEDP dataset of any form"""

    def __init__(self, data=None, axes=[], data_label='', log=[], attrs={} ): # Initialize object
        self.data = data
        self.axes = axes
        self.data_label =data_label
        self.log = log
        self.attrs = attrs
        
        if self.data is not None:
            #If data doesn't already have units, make it a dimensionless quantity
            if not isinstance(self.data, u.UnitBase ):
                self.data = self.data*u.Unit('' )
            # If axes haven't been created, make them up
            if len(self.axes) == 0:
                for i,n in enumerate(self.data.shape):
                    # Fill with a dimensionless unit
                    print('Creating axis')
                    a = {'name': 'ax' + str(i),   'axis' : np.arange(n)* u.Unit('' )}
                    self.axes.append(a)
            
            # Do some validation of the construction just in case
            if self.data.ndim != len(self.axes):
                print("ERROR: Number of axes does not match number of array dimensions! " + str(self.data.ndim) + ' != ' + str(len(self.axes))  )
                return None
            
            for i, ax in enumerate(self.axes):
                if len(ax['axis']) != self.data.shape[i]:
                    print("ERROR: Number of points in axis " + str( ax['name'] ) + ' does not match data shape: ' + str(self.data.shape)  )
                    return None

    def close(self):
        pass
    
    def copy(self):
        return copy.deepcopy(self)


    def readHDF(self, filepath):
        self.read_filepath = filepath
        with h5py.File(filepath, 'r') as f:   
            self.unpack(f)
        
    def saveHDF(self, filepath):
        self.save_filepath = filepath
        with h5py.File(filepath, 'w') as f:
            self.pack(f)


    def unpack(self, f):
        self.appendLog('Opened HDF file as NDF object: ' + self.read_filepath)

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
        
        self.appendLog('Saving NDF object as HDF ' + self.save_filepath)
        f['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    
    
    def getAxis(self, name):
        name = str(name)
        l = [ax['axis'] for ax in self.axes if ax['name'].lower().strip() == name.lower().strip()  ]
        if len(l) == 0:
            raise Exception("No axis found with name: " + str(name))
        elif len(l) > 1:
            raise Exception("Multiple axes found with name: " + str(name))
        else:
            return l[0]
    
    def getAxisInd(self, name):
        name = str(name)
        l = [i for i,ax in enumerate(self.axes) if ax['name'].lower().strip() == name.lower().strip()  ]
        if len(l) == 0:
            print("No axis found with name: " + str(name))
            return None
        elif len(l) > 1:
            print("Multiple axes found with name: " + str(name))
            return None
        else:
            return l[0]
    
    
    def avgDim(self, name ):
        ax_ind = self.getAxisInd(name)
        
        if ax_ind is None:
            print("DID NOT AVERAGE DIM: " + str(name))
            return None
        
        #Average the data
        self.data = np.average(self.data, axis=ax_ind)
        #Remove the axis from the the axes dictionary
        self.axes.pop(ax_ind)
        
        self.appendLog('Averaged over ' + name + ' axis')
        
        



    def collapseDim(self, name, value):
        # Find the axis index cooresponding to dim, and the axes it cooresponds to
        ax_ind = self.getAxisInd(name)
        
        if ax_ind is None:
            print("DID NOT COLLAPSE DIM: " + str(name))
            return None
        
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

        
        self.appendLog('Collapsed ' + name + ' axis to ' + name +  ' = ' + str(ax[ind]))
        
        
        
        
        
    def thinDim(self, name, bin=10):
        
        bin = int(bin)
        
        # Get the index that goes with this name
        ax_ind = self.getAxisInd(name)
        
        if ax_ind is None:
            print("DID NOT THIN DIM: " + str(name))
            return None

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
        
        if ax_ind is None:
            print("DID NOT CONVERT : " + str(name))
            return None
        
        ax = self.getAxis(name)
        
        old_unit = ax.unit    
        self.axes[ax_ind]['axis'] = ax.to(unit)

        self.appendLog('Converted ' + name + ' from ' + str(old_unit) + ' to ' + str( (self.axes[ax_ind]['axis']).unit ) )
        
        
        




    def plot(self, xrange=[None,None], yrange=[None,None], zrange=[None,None]):
        xname = self.axes[0]['name']
        xaxes = self.axes[0]['axis']
        
        if len(self.axes) == 1:
            print("Call simple 1D plotting routine")
            
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
            
            
            
            
    def __getattr__(self, key):
        key = key.lower().strip()
        axes_list = [ ax['name'] for ax in self.axes ]
        # If the axis is called, give the full astropy quantity object
        if key in axes_list:
            return self.getAxis(key)
        elif key in [ 'd' + s for s in axes_list] :
            # Strip off the 'd' to get the axis name
            k = key.replace('d', '')
            # Return the mean gradient as the step size
            return np.mean(  np.gradient( self.getAxis(k) ) )
       
        # Do some the same things for the dataset
        if key == 'data.unit':
            return self.data.unit
        elif key == 'data.value':
            return self.data.value
        
        #If you get this far, the key must be invalid
        print("Invalid Key: " + str(key))
        raise AttributeError
    
    
    #**************
    #HELPER METHODS
    #**************
    def appendLog(self, message):
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + str(message)
        self.log.append(entry)
        





class ndfxyz (ndf):
    """
    Class of gridded NDF datasets
    Inherits ndf
    """

    def unpack(self,f):
        ndf.unpack(self,f)
        
        self.x = self.getAxis('x')
        
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
    #sname = r"C:\Users\Peter\Desktop\TempData\testsave.h5"

    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    #fname = r"C:\Users\Peter\Desktop\TempData\run56_LAPD1_pos_raw.h5"
    #fname = r"C:\Users\Peter\Desktop\TempData\run102_PL11B_pos_raw.h5"
    
    
    
   # z = np.zeros([10,20])
   #a = [{'name':'first', 'axis': np.arange(10)*u.m},
   #       {'name':'sec', 'axis': np.arange(20)*u.mm}]
   # obj = ndf(data=z, axes=a)
    

    obj = ndf()
    obj.readHDF(fname)
    #obj.plot()
    #obj.data = obj.data*0 + 20 #Put in some fake data to test that it changes
    
    #print(obj.getAxis('reps'))
    
    obj.thinDim('t', bin=20)
    
    #obj2 = obj.copy()
    
    print(np.shape(obj.data))
    obj.avgDim('Reps')
    print(np.shape(obj.data))


    #print(obj.dtime.to(u.ns))
    
    #print(dir(obj))
    #print( obj.getAxis('time') )
    
    #print(np.shape(obj.data))
    #obj.collapseDim('x', 5*u.cm)
    
    #print(np.shape(obj.data))
    #obj.collapseDim('channels', 2)
    #print(np.shape(obj.data))
    
    #print(obj.x.unit)
    #obj.convertAxisUnit('x', u.mm)
    print(obj.dtime)
    
    #print(np.shape(obj.data))
    
    
    #obj.plot(xrange=[0,20 ])
    #obj.plot( xrange = [0, 40*u.us])
    
    obj.saveHDF(sname)


    obj.close()
    