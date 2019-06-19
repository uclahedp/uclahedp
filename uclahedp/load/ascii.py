# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:38:05 2019

@author: Peter
"""
import os
import numpy as np

import h5py

from uclahedp.tools import hdf as hdftools
from uclahedp.tools import csv as csvtools

def asciiToRaw(src, dest, delimiter=None, skip_header=None, ax=None, 
               axis_name=None, axis_unit=None, data_unit=None,
               run=None, probe=None, csv_dir=None):
     
     
     if csv_dir is not None:
          if run is not None or probe is not None:
               attrs = csvtools.getAllAttrs(csv_dir, run, probe)
 
     arr = np.genfromtxt(src.file, delimiter=delimiter, skip_header=skip_header )
     
     nelm, nchan = arr.shape

     
     if ax is None:
          axis = np.arange(nelm)
          axis_name = 'Indices'
     else:
          axis = arr[:,ax]
          axis_name = str(axis_name)
          
          outarr = np.zeros([nelm, nchan-1])
     
          i = 0
          for j in range(nchan):
               if j != ax:
                    outarr[:,i] = arr[:,j]
                    i +=1 
          arr = outarr
          

     
     #remove files if they already exist
     if os.path.exists(dest.file):
        os.remove(dest.file)
     
     #Create the destination file
     with h5py.File(dest.file, "a") as df:
          
        #Create the dest group, throw error if it exists
        if dest.group != '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
        grp = df[dest.group]
        
        #Initialize the output data array
        if 'data' in grp.keys():
            raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
            
        #Create the dataset + associated attributes
        grp.require_dataset("data", (nelm, nchan), np.float32, 
                            chunks=True, 
                            compression='gzip')
        grp['data'].attrs['unit'] = str(data_unit)
        
        grp['data'][:] = arr
        

        dimlabels = [str(axis_name), 'chan']
        
        grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
        
        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)

        #Create the axes
        grp.require_dataset(str(axis_name), (nelm,), np.float32, chunks=True )[:] = axis
        grp[str(axis_name)].attrs['unit'] = str(axis_unit)
        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ''
        
     return dest



if __name__ == "__main__":
     src = hdftools.hdfPath(os.path.join("F:","LAPD_Mar2018","INTERFEROMETER", "C1January 201700044.txt"))
     dest = hdftools.hdfPath(os.path.join("F:","LAPD_Mar2018","RAW", "44_Interferometer_Probe.hdf5"))
     
     csv_dir = os.path.join("F:","LAPD_Mar2018","METADATA", "CSV")
     run = 84
     probe = 'TonyInterferometer'
     
     asciiToRaw(src, dest, skip_header=5, ax=0, axis_name='time', 
                axis_unit='s', data_unit='V',
                csv_dir=csv_dir, run=run, probe=probe)