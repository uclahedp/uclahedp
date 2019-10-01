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


#Used for natural sorting filenames
import re

#Nice regex natural sorting algorithm found on stack overflow:
#https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def imgDirToRaw(run, probe, img_dir, dest, csv_dir, verbose=False):

     #Import attributes for this run/probe
     attrs = csvtools.getAllAttrs(csv_dir, run, probe)
  
         
     #Check for keys always required by this function
     req_keys = [ 'run_folder']
     csvtools.missingKeys(attrs, req_keys, fatal_error=True)
     
     run_folder = attrs['run_folder'][0]
     src = os.path.join(img_dir, run_folder)
     
     #Go through the directory and fine all the image files
     imgfiles = []
     for root, dirs, files in os.walk(src):  
         files = [f for f in files if f[0] != '.'] #Exclude files beginning in .
         for file in files:
             imgfiles.append(os.path.join(src, file))

     #Natural-sort the images by filename
     imgfiles = natural_sort(imgfiles)
     
     
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
                            chunks=(1, nxpx, nypx, 1), 
                            compression='gzip')
        grp['data'].attrs['unit'] = ''
        

        #Initialize time-remaining printout
        tr = util.timeRemaining(nframes, reportevery=5)


        #Actually put the images into the file
        for i,f in enumerate(imgfiles):
            tr.updateTimeRemaining(i)
            img = np.array(PIL.Image.open(f))
            
            img = np.reshape(img, [nxpx, nypx, nchan])
        
            #Rotate images
            for chan in range(nchan):
                img[:,:,chan] = np.rot90(img[:,:,chan], k=3)
            
            grp['data'][i, :,:,:] = img
         
            
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
    
     #data_dir = os.path.join("F:","LAPD_Jul2019")
     data_dir = os.path.join("/Volumes", "PVH_DATA","LAPD_Sept2019")
     
     
     
     img_dir = os.path.join(data_dir, 'PIMAX')
     dest = hdftools.hdfPath(os.path.join(data_dir ,"RAW", "run22.01_pimax4.hdf5"))
     csv_dir = os.path.join(data_dir,"METADATA")
     
     run = 22.01
     probe = 'pimax4'
     

     imgDirToRaw(run, probe, img_dir, dest, csv_dir, verbose=False)