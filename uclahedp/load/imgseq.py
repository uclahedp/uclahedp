# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:38:05 2019

@author: Peter
"""
import os
import numpy as np

import h5py, PIL

from uclahedp.tools import hdf as hdftools
from uclahedp.tools import csv as csvtools
from uclahedp.tools import util

def imgSeqToRaw(src, dest, run=None, probe=None, csv_dir=None, verbose=False):
     
     
     if csv_dir is not None and (run is not None or probe is not None):
         attrs = csvtools.getAllAttrs(csv_dir, run, probe)
     else:
         attrs = {}

     for root, dirs, files in os.walk(src):
         imgfiles = []
         imginds = []
         files = [f for f in files if f[0] != '.'] #Exclude files beginning in .
         for file in files:
             #We generate this integer for the purpose of sorting the images
             #This requires that the numbers in the filenames become more
             #specific as the name goes on, eg. run preceeds frame number
             strs= [s for s in file if s.isdigit()]
             sortint = int(''.join(strs))
             #Append both to the list of images
             imgfiles.append( os.path.join(src, file))
             imginds.append(sortint)
             
     #Sort the image list
     inds = np.argsort(np.array(imginds))
     imgfiles = [imgfiles[i] for i in inds]
     
     
     nframes = len(imgfiles)

  
     #remove files if they already exist
     if os.path.exists(dest.file):
        os.remove(dest.file)
     
     #Create the destination file
     with h5py.File(dest.file, "a") as df:
         
        #Assume all images are the same shape, load the first one to figure
        #out the array dimensions
        
        img = PIL.Image.open(imgfiles[0])
        nxpx, nypx = img.size
        #Bands will include the names of the different channels
        nchan = len(img.getbands())

        
        #Create the dest group, throw error if it exists
        if dest.group != '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
        grp = df[dest.group]
        
        #Initialize the output data array
        if 'data' in grp.keys():
            raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
            
        #Create the dataset + associated attributes
        grp.require_dataset("data", (nframes, nxpx, nypx, nchan), np.float32, 
                            chunks=True, 
                            compression='gzip')
        grp['data'].attrs['unit'] = ''
        

        #Initialize time-remaining printout
        tr = util.timeRemaining(nframes, reportevery=5)

        #Actually put the images into the file
        for i,f in enumerate(imgfiles):
            tr.updateTimeRemaining(i)
            img = PIL.Image.open(f)
            #TODO: FIX THIS SO IMAGES ARE ROTATED RIGHT!
            img = img.transpose(PIL.Image.TRANSPOSE) #Images are initially flipped
            grp['data'][i, :,:,:] = np.reshape(img, [nxpx, nypx, nchan])
         
            
        dimlabels = ['frames', 'xpixels', 'ypixels', 'chan']
        grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
        
        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)

        #Create the axes
        grp.require_dataset('frames', (nframes,), np.float32, chunks=True)[:] = np.arange(nframes)
        grp['frames'].attrs['unit'] = ''
        
        grp.require_dataset('xpixels', (nxpx,), np.float32, chunks=True)[:] = np.arange(nxpx)
        grp['xpixels'].attrs['unit'] = ''
        
        grp.require_dataset('ypixels', (nypx,), np.float32, chunks=True)[:] = np.arange(nypx)
        grp['ypixels'].attrs['unit'] = ''

        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ''
        
     return dest



if __name__ == "__main__":
     src = os.path.join('/Volumes','PVH_DATA','LAPD_Jul2019','PIMAX', 'run006')
     dest = hdftools.hdfPath(os.path.join('/Volumes','PVH_DATA','LAPD_Jul2019',"RAW", "run32_pimax006.hdf5"))
     
     #csv_dir = os.path.join("F:","LAPD_Jul2019","METADATA", "CSV")
     csv_dir = os.path.join("/Volumes", "PVH_DATA","LAPD_Jul2019","METADATA")
     
     run = 32
     probe = 'pimax006'
     
     imgSeqToRaw(src, dest, csv_dir=csv_dir, run=run, probe=probe)