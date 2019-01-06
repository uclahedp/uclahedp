#sdf class -> "standard data format"
#sdfarr -> labelled data array

import numpy as np
import h5py
import datetime
import matplotlib.pyplot as plt   
from astropy import units as u
import copy



class sdfarr:
    ##########################################################################
    # Basic Class Methods
    ##########################################################################
    def __init__(self, data=None, axes=[], data_label='', log=[]):
        """Initialize an sdfarr given keyword parameters"""
        self._setup(data, axes, data_label, log)
    
    def close(self):
        """Close this sdfarr object"""
        pass
        
    def copy(self):
        """Return a deep copy of this sdfarr object"""
        return copy.deepcopy(self)
    
    def readHDF(self, filepath, hdfpath = ''):
        """Read an sdfarr object from a location an an HDF file"""
        with h5py.File(filepath, 'r') as f: 
            g = f.require_group(hdfpath)
            self._unpack(g)
            try:
                self.appendLog('Opened HDF file: ' + filepath + ':' + hdfpath + ' as sdfarr' ) 
            except AttributeError:
                pass

    
    def saveHDF(self, filepath, hdfpath = ''):
        """Write an sdfarr object to an HDF file"""
        try:
            self.appendLog('Saving sdfarr in HDF file: ' + filepath + ':' + hdfpath)
        except AttributeError:
            pass
        with h5py.File(filepath) as f: 
            g = f.require_group(hdfpath)
            self._pack(g)
    
    
    
    def _setup(self, data, axes, data_label, log):
        """Setup the sdfarr, making some assumptions based on the keyword
        parameters given """
        self.data = data
        self.axes = axes
        self.data_label = data_label
        self.log = log

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
                    
                    
                    
    def _unpack(self, g):
        """Unpack the object from an HDF file group g"""
        self.data =  g['data'][:] *  u.Unit( g['data'].attrs['unit'], parse_strict = 'warn' )
        self.data_label = g['data'].attrs['label']
        self.log = g['log'][:]
        self.axes = []
        for i in range(self.data.ndim):
            name = 'ax' + str(i)
            label = g[name].attrs['label']
            unit = u.Unit(g[name].attrs['unit'], parse_strict = 'warn' )
            ax = g[ name ][:] * unit 
            a = {'label':label, 'axis':ax}
            self.axes.append(a)
            
    def _pack(self, g):
        """Pack the object into an HDF file group g"""
        #Delete any keys currently existing in the dataset path
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
            
        g['log'] = [s.encode('utf-8') for s in self.log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
    def appendLog(self, message):
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + str(message)
        self.log.append(entry)
        
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
        
    ##########################################################################
    # Array operations
    ##########################################################################
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
        return self.axes[self.getAxisInd(label)]['axis']

    def avgDim(self, label ):
        """Average over the axis specified by the given label"""
        ax_ind = self.getAxisInd(label)
        #Average the data
        self.data = np.average(self.data, axis=ax_ind)
        #Remove the axis from the the axes dictionary
        self.axes.pop(ax_ind)
        self.appendLog('Averaged over ' + label + ' axis')
        
    def collapseDim(self, label, value):
        """Collapse the dimension specified to the element closest to value
        For example, collapseDim('x', 5) collapses the x dimension to the value 5
        """
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
        self.appendLog('Collapsed ' + label + ' axis to ' + label +  ' = ' + str(ax[ind]))
        
        
    def thinDim(self, label, bin=10):
        #Coerrce bin to type int
        bin = int(bin)
        # Get the index that goes with this name
        ax_ind = self.getAxisInd(label)
        ax = self.getAxis(label)
        # Get the current shape vector
        shape = list( self.data.shape )
        # Create a shape array for a single slice
        slice_shape = list( self.data.shape )
        slice_shape[ax_ind] = 1
      
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
            # Reshape here makes sure the dimensions match
            arr[ tuple(pslice) ] = np.reshape(s, slice_shape)
            
            #Now do the same thing for the cooresponding axis
            s = np.take(ax, indrange)
            s = np.average(s)
            np.put(newax, i, s)
        # Put the data back into the array
        self.data = arr*self.data.unit
        # Put the new axis back too
        self.axes[ax_ind]['axis'] = newax * self.axes[ax_ind]['axis'].unit
        
        self.appendLog("Thinned " + label + " axis with bin size " + str(bin))
        
    def convertAxisUnit(self, label, unit):
        unit = self._coerceAstropyUnit(unit)
    
        ax_ind = self.getAxisInd(label)
        ax = self.getAxis(label)

        old_unit = ax.unit    
        self.axes[ax_ind]['axis'] = ax.to(unit)
        
        self.appendLog('Converted ' + label + ' axis units from ' + str(old_unit) + ' to ' + str( (self.axes[ax_ind]['axis']).unit ) )
        
    def convertDataUnit(self, unit):
        unit = self._coerceAstropyUnit(unit)
        old_unit = self.data.unit    
        self.data = self.data.to(unit)
        self.appendLog('Converted data units from ' + str(old_unit) + ' to ' + str( self.data.unit )   )
    
    ##########################################################################
    # Very Basic Plotter Methods
    ##########################################################################
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
    
    
    ##########################################################################
    # Helper Methods
    ##########################################################################
    def _coerceAstropyUnit(self, unit):
        """Convert a string to an astropy unit if possible
        If it is already an astropy unit, just return it.
        If the conversion is impossible, raise an exception.
        """
        if isinstance(unit, str ):
            try:
                unit = u.Unit( unit )
            except ValueError:
                raise("Reqested axis unit is not recognized by astropy.units: " + str(unit) )

        if not isinstance(unit, u.UnitBase):
            raise("Reqested axis unit is not recognized by astropy.units: " + str(unit) )
            
        return unit
    
    
        
        
        
        
        
if __name__ == "__main__":
    sname = r"C:\Users\Peter\Desktop\TempData\testsave.h5"
    
    fname = r"C:\Users\Peter\Desktop\TempData\run56_LAPD1_pos_raw.h5"
    
    z = np.ones([10,15,200,3])*u.V
    a = [{'label':'x', 'axis': np.arange(10)*u.cm},
          {'label':'y', 'axis': np.arange(15)*u.mm},
          {'label':'z', 'axis': np.arange(200)*u.mm},
          {'label':'rep', 'axis': np.arange(3)}]
   
    attrs = {'one': 1, 'two' : [1,2,3] }
    
    obj = sdfarr(data=z, axes=a)
    
    obj.convertAxisUnit('z', u.mm)
    obj.convertDataUnit('mV')

    
    obj.avgDim('rep')
    obj.collapseDim('x', 3*u.cm)
    obj.thinDim('z', bin=5)
    
    obj.saveHDF(sname, '/test/dir/')
    obj.close()