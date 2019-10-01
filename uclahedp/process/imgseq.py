#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:50:31 2019

@author: peter
"""

import os, h5py
import numpy as np


from uclahedp.tools import hdf as hdftools
from uclahedp.tools import csv as csvtools
from uclahedp.tools import util

import astropy.units as u



def imgSeqRawToFull(src, dest):
    with h5py.File(src.file, 'r') as sf:
             
            #Get the datagroup
            srcgrp = sf[src.group]
            
            #Create dictionary of attributes
            attrs = hdftools.readAttrs(srcgrp)
            
            #Check for keys always required by this function
            req_keys = [ 'dt']
            
            csvtools.missingKeys(attrs, req_keys, fatal_error=True)
            
            nframes, nxpx, nypx, nchan = srcgrp['data'].shape
            
            
            #Convert dt
            dt = (attrs['dt'][0]*u.Unit(attrs['dt'][1])).to(u.s).value
            
            
            
            #Reps is assumed to be 1 unless otherwise set
            if 'nreps' in attrs.keys():
                nreps = attrs['nreps'][0]
            else:
                nreps = 1
                
                
            nti = int(nframes/nreps)
            
            #t0 is the time of the first frame in the set
            if 't0' in attrs.keys():
                t0 = (attrs['t0'][0]*u.Unit(attrs['t0'][1])).to(u.s).value
            else:
                t0 = 0
            
            #Laser t0 is the time when the laser fires
            #Time array will be shifted so this time is zero
            if 'laser_t0' in attrs.keys():
                laser_t0 = (attrs['laser_t0'][0]*u.Unit(attrs['laser_t0'][1])).to(u.s).value
            else:
                laser_t0 = 0
                
                
            #dxdp is the pixel spacing in cm/px
            if'dxdp' in attrs.keys():
                dxdp = (attrs['dxdp'][0]*u.Unit(attrs['dxdp'][1])).to(u.cm).value
            else:
                dxdp = None
                
            if'dydp' in attrs.keys():
                dydp = (attrs['dydp'][0]*u.Unit(attrs['dydp'][1])).to(u.cm).value
            else:
                dydp = None
                
            
            if'x0' in attrs.keys():
                x0 = (attrs['x0'][0]*u.Unit(attrs['x0'][1])).to(u.cm).value
            else:
                x0 = 0
                
            if'y0' in attrs.keys():
                y0 = (attrs['y0'][0]*u.Unit(attrs['y0'][1])).to(u.cm).value
            else:
                y0 = 0
                
            
            with h5py.File(dest.file, 'a') as df:
                destgrp = df.require_group(dest.group)
                
                destgrp.require_dataset("data", (nti, nxpx, nypx, nreps, nchan), np.float32, 
                            chunks=(1, nxpx, nypx, 1, 1), 
                            compression='gzip')
                destgrp['data'].attrs['unit'] = ''
                
                #Initialize time-remaining printout
                tr = util.timeRemaining(nti, reportevery=5)

                #Actually put the images into the file
                for i in range(nti):
                    tr.updateTimeRemaining(i)
                    
                    a = i*nreps
                    b = (i+1)*nreps
                    
                    #print(str(a) + ":" + str(b))
                    
                    #Copy, re-shape, and write data to array
                    arr = srcgrp['data'][a:b,...]

                    
                    arr = np.moveaxis(arr, 0, 2)

                    
                    #arr = np.reshape(arr, [nreps, nxpx, nypx, nchan])
                    destgrp['data'][i,...] = arr
                
                
                
                
                #Write the attrs dictioanry into attributes of the new data group
                hdftools.writeAttrs(attrs, destgrp)
                
            
                dimlabels = []
                
             
                time = np.arange(nti)*dt + t0 - laser_t0
                destgrp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = time
                destgrp['time'].attrs['unit'] = 's'
                dimlabels.append('time')
                
                
                if dxdp is not None:
                    xaxis = np.arrange(nxpx)*dxdp + x0
                    destgrp.require_dataset('xaxis', (nxpx,), np.float32, chunks=True)[:] = xaxis
                    destgrp['xaxis'].attrs['unit'] = 'cm'
                    dimlabels.append('xaxis')
                else:
                    destgrp.require_dataset('xpixels', (nxpx,), np.float32, chunks=True)[:] = np.arange(nxpx)
                    destgrp['xpixels'].attrs['unit'] = ''
                    dimlabels.append('xpixels')
        
        
                if dydp is not None:
                    yaxis = np.arrange(nypx)*dydp + y0
                    destgrp.require_dataset('yaxis', (nypx,), np.float32, chunks=True)[:] = yaxis
                    destgrp['yaxis'].attrs['unit'] = 'cm'
                    dimlabels.append('yaxis')
                else:
                    destgrp.require_dataset('ypixels', (nypx,), np.float32, chunks=True)[:] = np.arange(nypx)
                    destgrp['ypixels'].attrs['unit'] = ''
                    dimlabels.append('ypixels')
                    
                destgrp.require_dataset('reps', (nreps,), np.float32, chunks=True)[:] = np.arange(nreps)
                destgrp['reps'].attrs['unit'] = ''
                dimlabels.append('reps')
                    
                destgrp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
                destgrp['chan'].attrs['unit'] = ''
                dimlabels.append('chan')
                
                destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
              
            

                
                
if __name__ == "__main__":
    
     #data_dir = os.path.join("F:","LAPD_Jul2019")
     data_dir = os.path.join("/Volumes", "PVH_DATA","LAPD_Sept2019")
     
     src = hdftools.hdfPath(os.path.join(data_dir, 'RAW', 'run15.01_pimax4.hdf5'))
     dest = hdftools.hdfPath(os.path.join(data_dir ,"FULL", "run15.01_pimax4_full.hdf5"))

     imgSeqRawToFull(src, dest)