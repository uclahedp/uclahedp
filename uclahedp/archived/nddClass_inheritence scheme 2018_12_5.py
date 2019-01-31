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

ndd -> "N-Dimensional Dataset
"""


class ndd_base:
    def __init__(self):
        pass
    
    def close(self):
        pass
    
    def copy(self):
        return copy.deepcopy(self)
    
    def pack(self, g):
        pass
    
    def unpack(self, g):
        pass


    def readHDF(self, filepath, hdfpath = ''):
        with h5py.File(filepath, 'r') as f: 
            g = f.require_group(hdfpath)
            self.unpack(g)
            
            try:
                self.log.append('Opened HDF file: ' + filepath + ':' + hdfpath + ' as ndd object ' + str(type(self)) ) 
            except AttributeError:
                pass

    
    def saveHDF(self, filepath, hdfpath = ''):
        try:
            self.log.append('Saving ndd object ' + str(type(self)) + ' in HDF file: ' + filepath + ':' + hdfpath)
        except AttributeError:
            pass
        
        with h5py.File(filepath) as f: 
            g = f.require_group(hdfpath)
            
            self.pack(g)
            


class ndd_log(ndd_base):
    def __init__(self, log=[] ):
        self.setup(log)
        
    def setup(self, log):
        self.log = log
        
    def append(self, message):
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + str(message)
        self.log.append(entry)
        
    def getLog(self):
        return self.log
    
    def pack(self, g):
        
        if 'log' in g.keys():
            del(g['log'])

        g['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
         
        
    def unpack(self, g):
        self.log = g['log'][:]
         
        

        
        
        
class ndd_attrs (ndd_base):
    def __init__(self, attrs={}, log=None ):
        self.setup(attrs)
        
        # Log is an optional argument
        # if a log object is used, that will be associated with this object
        # if a list is used, that will be used to create a log object inside this object
        # if nothing is given, no log will be added.
        if isinstance(log, ndd_log):
            self.log = log
        else:
            self.log = ndd_log(log)

    
    def setup(self, attrs):
        self.attrs = attrs
        
        
    def __getattr__(self, key):
        if key == 'attrs':
            return self.attrs


    def unpack(self, f):
        self.attrs = {}
        keys = list(f.attrs)
        for key in keys:
            self.attrs[key] = f.attrs[key]
            
            
    def pack(self, f):
        for key in self.attrs.keys():
            f.attrs[key] =self.attrs[key]

        
        
        
        

class ndd_arr (ndd_base):
    def __init__(self, data=None, axes=[], data_label='', log=None): # Initialize object
        self.setup(data, axes, data_label)
        
        # Log is an optional argument
        # if a log object is used, that will be associated with this object
        # if a list is used, that will be used to create a log object inside this object
        # if nothing is given, no log will be added.
        
        if isinstance(log, ndd_log):
            self.log = log
        else:
            self.log = ndd_log(log)
            

        
        
    def setup(self, data, axes, data_label):
        self.data = data
        self.axes = axes
        self.data_label = data_label
        
        
        if self.data is not None:
            #If data doesn't already have units, make it a dimensionless quantity
            if not isinstance(self.data, u.UnitBase ):
                self.data = self.data*u.Unit('' )
            # If axes haven't been created, make them up
            if len(self.axes) == 0:
                for i,n in enumerate(self.data.shape):
                    # Fill with a dimensionless unit
                    print('Creating axis')
                    a = {'label': 'ax' + str(i),   'axis' : np.arange(n)* u.Unit('' )}
                    self.axes.append(a)
            
            # Do some validation of the construction just in case
            if self.data.ndim != len(self.axes):
                raise("ERROR: Number of axes does not match number of array dimensions! " + str(self.data.ndim) + ' != ' + str(len(self.axes))  )

            for i, ax in enumerate(self.axes):
                if len(ax['axis']) != self.data.shape[i]:
                    raise("ERROR: Number of points in axis " + str( ax['label'] ) + ' does not match data shape: ' + str(self.data.shape)  )
                    




    def unpack(self, g):
        

        self.data =  g['data'][:] *  u.Unit( g['data'].attrs['unit'], parse_strict = 'warn' )
        self.data_label = g['data'].attrs['label']

        self.axes = []
        for i in range(self.data.ndim):
            name = 'ax' + str(i)
            label = g[name].attrs['label']
            unit = u.Unit(g[name].attrs['unit'], parse_strict = 'warn' )
            ax = g[ name ][:] * unit 
            a = {'label':label, 'axis':ax}
            self.axes.append(a)


        
    def pack(self, g):
        
        for key in g.keys():
            del(g[key])

        g['data'] = self.data
        g['data'].attrs['label'] = self.data_label
        g['data'].attrs['unit'] = str(self.data.unit)
    
        for i, ax in enumerate(self.axes):
            name = 'ax' + str(i)
            g[name] = ax['axis'].value
            g[name].attrs['label'] = ax['label']
            g[name].attrs['unit'] = str( ax['axis'].unit )

    

    
    def getAxisInd(self, label):
        """Return the index of the axis cooresponding to the given axis label"""
        label = str(label)
        l = [i for i,ax in enumerate(self.axes) if ax['label'].lower().strip() == label.lower().strip()  ]
        if len(l) == 0:
            raise("No axis found with label: " + str(label))
        elif len(l) > 1:
           raise("Multiple axes found with label: " + str(label))
        else:
            return l[0]
        
     def getAxis(self, label):
        """ Return the axis array cooresponding to a particular axis label """ 
        return self.axes[self.getAxisInd(label)]
    
    
    def avgDim(self, label ):
        ax_ind = self.getAxisInd(label)
        

        #Average the data
        self.data = np.average(self.data, axis=ax_ind)
        #Remove the axis from the the axes dictionary
        self.axes.pop(ax_ind)
     
        self.log.append('Averaged over ' + label + ' axis')
        
        



    def collapseDim(self, label, value):
        # Find the axis index cooresponding to dim, and the axes it cooresponds to
        ax_ind = self.getAxisInd(label)
      
        ax = self.getAxis(label)

        #If value isn't already an astropy quantity, assume it has the same units as the axis.
        if not isinstance(value, u.Quantity ):
            value = value * ax.unit

        #Find the index along the axis that is closes to 'value'
        ind = np.where(abs(ax - value) == abs(ax - value).min())[0][0]

        #Collapse the dimension in the array
        self.data = np.take(self.data, ind, axis=ax_ind)
        #Delete the axis from the object
        self.axes.pop(ax_ind)

        
        self.log.append('Collapsed ' + label + ' axis to ' + label +  ' = ' + str(ax[ind]))
        
        
        
        
        
    def thinDim(self, label, bin=10):
        
        bin = int(bin)
        
        # Get the index that goes with this name
        ax_ind = self.getAxisInd(label)
        ax = self.getAxis(label)
        
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
        
        

    def convertAxisUnit(self, label, unit):
        if isinstance(unit, str ):
            try:
                unit = u.Unit( unit )
            except ValueError:
                raise("Reqested axis unit is not recognized by astropy.units: " + str(unit) )

        if not isinstance(unit, u.UnitBase):
            raise("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
        
        ax_ind = self.getAxisInd(label)
        ax = self.getAxis(label)
        
        old_unit = ax.unit    
        self.axes[ax_ind]['axis'] = ax.to(unit)

        self.log.append('Converted ' + label + ' from ' + str(old_unit) + ' to ' + str( (self.axes[ax_ind]['axis']).unit ) )
        
        
        




    def plot(self, xrange=[None,None], yrange=[None,None], zrange=[None,None]):
        xname = self.axes[0]['label']
        xaxes = self.axes[0]['axis']
        
        if len(self.axes) == 1:
            print("Call simple 1D plotting routine")
            
            #Make plot
            plt.plot(xaxes, self.data)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(xname + ' (' + str(xaxes.unit) + ')'  )
            plt.ylabel(self.data_label.title() + ' (' + str(self.data.unit) + ')'  )
            plt.show()
            plt.close()
            
        elif len(self.axes) == 2:
            print("Call simple 2D contour plotting routine")
            yname = self.axes[1]['label']
            yaxes = self.axes[1]['axis']
            

            plt.figure()
            plt.contourf(xaxes, yaxes, self.data.T)
            plt.axis([xrange[0], xrange[1] , yrange[0],  yrange[1]])
            plt.xlabel(xname + ' (' + str(xaxes.unit) + ')'  )
            plt.ylabel(yname + ' (' + str(yaxes.unit) + ')'  )
            plt.title(self.data_label.title() + ' (' + str(self.data.unit) + ')'  )
            plt.show()
            plt.close()
            
        elif len(self.axes) == 3:
            print("Ugh call a 3D plotting routine yuck")
        
        else:
            print("STOP TRYING TO VISUALIZE 4+ SPATIAL DIMENSIONS")
            
            
            
            
    def __getattr__(self, key):
        key = key.lower().strip()
        axes_list = [ ax['label'] for ax in self.axes ]
        # If the axis is called, give the full astropy quantity object
        if key in axes_list:
            return self.getAxis(key)
        elif key in [ 'd' + s for s in axes_list] :
            # Strip off the 'd' to get the axis name
            k = key.replace('d', '')
            # Return the mean gradient as the step size
            return np.mean(  np.gradient( self.getAxis(k) ) )
       

        
        #If you get this far, the key must be invalid
        print("Invalid Key: " + str(key))
        raise AttributeError
    

        
        
class ndd(ndd_base):
    """
    ndd is an ndd_arr that has an attached dictionary of attributes
    """
    def __init__(self, data=None, axes=[], data_label='', log=[], attrs={} ): # Initialize object
        self.log = ndd_log( log)
        self.attrs = ndd_attrs( attrs, log=self.log)
        self.arr = ndd_arr(data, axes, data_label, log=self.log)

    def unpack(self, f):
        self.log.unpack(f) #First so all messages get recorded
        self.attrs.unpack(f)
        self.arr.unpack(f)
        pass
        
        
         
    def pack(self, f):
        self.attrs.pack(f)
        self.arr.pack(f)
        self.log.pack(f) #Last so all messages get recorded
        pass

  
    
    def __getattr__(self, key):
        axes_list = [ ax['label'] for ax in self.arr.axes ]
        if key == "data":
            return self.arr.data
        elif key in axes_list:
            return self.arr.getAxis(key)
        elif key in [ 'd' + s for s in axes_list] :
            # Strip off the 'd' to get the axis name
            k = key.replace('d', '')
            # Return the mean gradient as the step size
            return np.mean(  np.gradient( self.arr.getAxis(k) ) )
        elif key == "ax_labels":
            return [ax['label'] for ax in self.arr.axes]
        elif key == "attr_keys":
            return [k for k in self.attrs.attrs.keys()]

        
    
    # Allow some ndd_arr methods to be directly callable on the ndd object
    # is this a good idea?
    def avgDim(self, label):
        self.arr.avgDim(label)
    def collapseDim(self, label, value):
        self.arr.collapseDim(label, value)
    def thinDim(self, label, bin):
        self.arr.thinDim(label, bin)
    def plot(self, xrange, yrange,zrange):
        self.arr.plot(xrange, yrange, zrange)
 
        

        
if __name__ == "__main__":
    sname = r"/Users/peter/Desktop/tempdata/testsave2.h5"
    #sname = r"C:\Users\Peter\Desktop\TempData\testsave.h5"

    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_full.h5"
    fname = r"/Users/peter/Desktop/tempdata/run56_LAPD1_pos_raw.h5"
    #fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_full.h5"
    
    #fname = r"C:\Users\Peter\Desktop\TempData\run56_LAPD1_pos_raw.h5"
    #fname = r"C:\Users\Peter\Desktop\TempData\run102_PL11B_pos_raw.h5"
    


    z = np.ones([10,20])*u.V
    a = [{'label':'x', 'axis': np.arange(10)*u.m},
          {'label':'y', 'axis': np.arange(20)*u.mm} ]
   
    attrs = {'one': 1, 'two' : [1,2,3] }
    a2 = {'another': 1}

    obj = ndd(data=z, axes=a, attrs=attrs)
    obj.log.append("message1")
    obj.saveHDF(sname, '/run56/LAPD1/')
    obj.close()
    
    
    obj = ndd(data=z*2, axes=a, attrs=attrs)
    obj.log.append("message2")
    obj.saveHDF(sname, '/run56/LAPD10/')
    obj.close()
    
    
    obj = ndd(data=z*3, axes=a, attrs=attrs)
    obj.log.append("Overwriting with 3's")
    obj.avgDim('x')
    obj.saveHDF(sname, '/run56/LAPD1/')
    obj.close()
    
    
    
    obj = ndd_attrs(attrs = {'test add to atters':1})
    obj.saveHDF(sname, '/run56/LAPD1/')
    obj.close()
    
    obj = ndd_attrs(a2)
    obj.saveHDF(sname, '/run56/')
    obj.close()
    
    obj = ndd()
    obj.readHDF(sname,  '/run56/LAPD1/' )
    print(obj.attr_keys)
    obj.close()
    
    
    
    
    
    
    
    

    