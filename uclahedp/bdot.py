#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bdot.py: BDOT analysis package

Created on Wed Nov 28 13:37:21 2018

@author: peter
"""

import csvtools
import hdftools
import util
import tdiode
import postools

import numpy as np
import os
import h5py
from scipy.signal import detrend as detrend
import astropy.units as u
import time




def bdot_raw_to_full(src, dest, tdiode_hdf=None, grid=False, verbose=False):
    """ Integrates bdot data, calibrates output using information about the probe.
        Corrects for probe angle based on which drive is being used.

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing bdot data
            
        dest: hdfPath object
            Path string to location processed bdot data should be written out

        tdiode:  hdfPath object
            Path to a tdiode full hdf5 dataset containing diode t0 times and
            an array of boolean 'bad shot' flags.
            
        grid: Boolean
            If grid is true, output will be written in cartesian grid array
            format, eg. [nti, nx, ny, nz, nreps, nchan]. Otherwise, output will
            be in [nshots, nti, nchan] format


    Returns
    -------
       None
    """ 
    # ******
    # Load data from the raw HDF file
    # ******
    with h5py.File(src.file, 'r') as sf:
        
        kdict = {}
        srcgrp = sf[src.group]
        
        #src_run_grp = srcgrp.parent
  
        
        #Check for keys always required by this function
        req_keys = ['brd1','brd2','brd3', 'chan1','chan2', 'chan3', 
                    'xarea', 'yarea', 'zarea',
                    'xatten', 'yatten', 'zatten', 'gain',
                    'xpol', 'ypol', 'zpol', 'roll', 
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z',
                    'dt', 'nturns']
       
        motion = None
        if  'pos' in srcgrp:
            pos = srcgrp['pos'][:] #Read the entire array in
            req_keys = req_keys + ['motion_format']
            #If pos array exists, there are keywords required for that too.
            if srcgrp.attrs['motion_format'] == 'fixed_rotation':
                req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
            else:
                raise ValueError("Motion format unrecognized: " + str(srcgrp['pos'].attrs['motion_format']) )
            motion = srcgrp.attrs['motion_format']
        else:
            #If no position information is given, a single explicit position
            #is required. 
            req_keys = req_keys + ['xpos', 'ypos', 'zpos']
            grid = False #Can't grid data if there's no pos array!
            
            
        #Process the required keys, throwing an error if any cannot be found
        missing_keys = []
        for k in req_keys:
            if k in srcgrp.attrs:
                kdict[k] = ( csvtools.fixType(srcgrp.attrs[k][0]), srcgrp.attrs[k][1])
            else:
                missing_keys.append(k)
        if len(missing_keys) > 0:
            raise ValueError("Missing required keys! ->" + str(missing_keys))
        

        nshots, nti, nchan = srcgrp['data'].shape
        
        
        #If keywords specify gridded output, run postools.grid 
        if grid:
            shotgridind, xaxes, yaxes, zaxes = postools.gridShotIndList(pos, precision=.1)
            nx,ny,nz = len(xaxes), len(yaxes), len(zaxes)
            nreps =int( np.floor(nshots/(nx*ny*nz)))
            

        if tdiode_hdf is not None:
            if verbose:
                print("Loading tdiode array from file.")
            #Get an array of all the t0 indices
            t0indarr = tdiode.calcT0ind(tdiode_hdf)
            #Get an array of all the good shots and badshots (indices)
            badshots, goodshots = tdiode.findBadShots(tdiode_hdf)
            #Replace any bad shots with a standin value for now
            t0indarr[badshots] = int(np.median(t0indarr[goodshots]))
            
            min_t0ind = np.min(t0indarr[goodshots])
            max_t0shift = np.max(t0indarr[goodshots]) - min_t0ind
            nti = nti - max_t0shift
            
            
            

        with h5py.File(dest.file, 'a') as df:
            
            #Throw an error if this group already exists
            if dest.group is not '/' and dest.group in df.keys():
                raise hdftools.hdfGroupExists(dest)
            
            destgrp = df.require_group(dest.group)

            
            #Copy over attributes
            """
            for k in src_run_grp.attrs.keys():
                dest_run_grp.attrs[k] = src_run_grp.attrs[k]
            for k in srcgrp.attrs.keys():
                destgrp.attrs[k] = srcgrp.attrs[k]
            """
            for k in srcgrp.attrs.keys():
                destgrp.attrs[k] = srcgrp.attrs[k]
            
            #Throw an error if this dataset already exists
            if 'data' in destgrp.keys():
                    raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
                    
            #Create the dataset 'data' appropriate to whether or not output
            #data will be gridded
            if verbose:
                print("Creating 'data' group in destination file")
            if grid:
                destgrp.require_dataset('data', (nti, nx, ny, nz, nreps, nchan), np.float32, chunks=(20000,1,1,1,1,1), compression='gzip')
            else:
                destgrp.require_dataset('data', (nshots, nti, nchan), np.float32, chunks=(1, 20000, 1), compression='gzip')
            
            #Load the time vector
            t = srcgrp['time']
            #If a timing diode is being applied, correct the time vector here.
            if tdiode_hdf is not None:
                t = t[0:nti] - t[min_t0ind]

            atten = np.array([kdict['xatten'][0],kdict['yatten'][0],kdict['zatten'][0]])
            # Create an array of calibration coefficents
            if kdict['xatten'][1] == 'dB':
                print("Converting dB to x")
                atten = np.power([10,10,10], atten/20.0) # Convert from decibels
           
            # dt -> s
            dt = ( kdict['dt'][0]*u.Unit(kdict['dt'][1])).to(u.s).value

            # area : mm^2 -> m^2
            xarea = (kdict['xarea'][0]*u.Unit(kdict['xarea'][1])).to(u.m ** 2).value
            yarea = (kdict['yarea'][0]*u.Unit(kdict['yarea'][1])).to(u.m ** 2).value
            zarea = (kdict['zarea'][0]*u.Unit(kdict['zarea'][1])).to(u.m ** 2).value

            gain = kdict['gain'][0]
            nturns = kdict['nturns'][0]

            xcal = 1.0e4*dt*atten[0]/gain/(nturns*xarea)
            ycal = 1.0e4*dt*atten[1]/gain/(nturns*yarea)
            zcal = 1.0e4*dt*atten[2]/gain/(nturns*zarea)
            
            xpol = kdict['xpol'][0]
            ypol = kdict['ypol'][0]
            zpol = kdict['zpol'][0]
            
            
            #Initialize some variables to use in time-remaining printout
            tstart  = time.time()
            tperstep = []
            nstepsreport = int(nshots/100.0)
            
            if verbose:
                print("Beginning processing data shot-by-shot.")
            
            #Chunking data processing loop limits memory usage
            for i in range(nshots):
                #Update the time-per-shot time estimator
                #Over time this should give more accurate run-time
                #predictions
                nowtime = time.time()
                tperstep.append(nowtime-tstart)
                tstart = nowtime
                    
                if verbose and i % nstepsreport == 0 and i > 0:
                    tremain = np.mean(tperstep)*(nshots - i)
                    print(str(i) + '/' + str(nshots) + ' complete, ' + 
                          util.timeFormat(tremain) + ' remaining' )
                
                

                #If a tdiode hdf was supplied, apply the correction here.
                if tdiode_hdf is not None:
                    ta = t0indarr[i] - min_t0ind
                    tb = ta + nti
                   # print('t0ind:' + str(t0indarr[i]) +  ' ta: ' + str(ta) + 
                   #    ', tb: ' + str(tb) + ', nti: ' + str(nti) + 
                   #    ', tb-ta: ' +  str(tb-ta))
                else:
                    ta = 0
                    tb = -1
                
                    
                bx = srcgrp['data'][i,ta:tb, 0]
                by = srcgrp['data'][i,ta:tb, 1]
                bz = srcgrp['data'][i,ta:tb, 2]
                
                bx = np.cumsum( detrend(bx) )*xcal*xpol
                by = np.cumsum( detrend(by) )*ycal*ypol
                bz = np.cumsum( detrend(bz) )*zcal*zpol
                
     
                if motion == 'fixed_rotation':
                    x,y,z = srcgrp['pos'][i, :]
                    rx, ry, rz = kdict['rot_center_x'][0],kdict['rot_center_y'][0],kdict['rot_center_z'][0]
                    pitch = np.arctan( (y-ry) / (x-rx) ) 
                    yaw = np.arctan( (z-rz) / (x-rx) ) 
                    
                    roll, unit = kdict['roll']
                    if unit != 'rad':
                        np.radians(roll)
              
                    #Matrix is the first Tait-Bryan matrix XZY from https://en.wikipedia.org/wiki/Euler_angles
                    #1 -> roll
                    #2 -> pitch
                    #3 -> yaw
                    bx = (np.cos(pitch)*np.cos(yaw)*bx - 
                        np.sin(pitch)*by  + 
                        np.cos(pitch)*np.sin(yaw)*bz)
                    
                    by =  ((np.sin(roll)*np.sin(yaw) + np.cos(roll)*np.cos(yaw)*np.sin(pitch))*bx +
                           np.cos(roll)*np.cos(pitch)*by  +
                           (np.cos(roll)*np.sin(pitch)*np.sin(yaw) - np.cos(yaw)*np.sin(roll))*bz)
                    
                    bz =  ((np.cos(yaw)*np.sin(roll)*np.sin(pitch) - np.cos(roll)*np.sin(yaw))*bx + 
                           np.cos(pitch)*np.sin(roll)*by  +
                           (np.cos(roll)*np.cos(yaw) + np.sin(roll)*np.sin(pitch)*np.sin(yaw))*bz)
                    
                if grid:

                    xi = shotgridind[i, 0]
                    yi = shotgridind[i, 1]
                    zi = shotgridind[i, 2]
                    repi = shotgridind[i, 3]

                    destgrp['data'][:, xi, yi, zi, repi, 0] = bx
                    destgrp['data'][:, xi, yi, zi, repi, 1] = by
                    destgrp['data'][:, xi, yi, zi, repi, 2] = bz
                else:
                    destgrp['data'][i,:, 0] = bx
                    destgrp['data'][i,:, 1] = by 
                    destgrp['data'][i,:, 2] = bz                      
            
            
            
            
            if verbose:
                print("Writing axes to destination file")
            
            if motion is not None:
                #Add the other axes and things we'd like in this file
                destgrp.require_dataset('pos', (nshots, 3), np.float32, chunks=True)[:] = srcgrp['pos'][:]
                for k in srcgrp['pos'].attrs.keys():
                    destgrp['pos'].attrs[k] = srcgrp['pos'].attrs[k]

            if grid:
                destgrp.require_dataset('xaxes', (nx,), np.float32, chunks=True)[:] = xaxes
                destgrp['xaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']
                
                destgrp.require_dataset('yaxes', (ny,), np.float32, chunks=True)[:] = yaxes
                destgrp['yaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']
                
                destgrp.require_dataset('zaxes', (nz,), np.float32, chunks=True)[:] = zaxes
                destgrp['zaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']
                
                dimlabels = ['time', 'xaxes', 'yaxes', 'zaxes', 'reps', 'chan']
                destgrp['data'].attrs['shape'] = [nti, nx, ny, nz, nreps, nchan]
            else:
                destgrp.require_dataset('shots', (nshots,), np.int32, chunks=True)[:] = srcgrp['shots'][:]
                destgrp['shots'].attrs['unit'] = srcgrp['shots'].attrs['unit']
                dimlabels = ['shots', 'time', 'chan']
                destgrp['data'].attrs['shape'] = [nshots, nti, nchan]
            
            
            
            destgrp.require_dataset('chan', (nchan,), np.int32, chunks=True)[:] = srcgrp['chan'][:]
            destgrp['chan'].attrs['unit'] = srcgrp['chan'].attrs['unit']

            
            try:
                destgrp.require_dataset('time', (nti,), np.float32, chunks=True)
            except ValueError:
                destgrp['time'].resize( (nti,))
                #destgrp.require_dataset('time', (nti,), np.float32, chunks=(1, nti, 1) )
                
                
        
            destgrp['time'][:] = t
            destgrp['time'].attrs['unit'] = srcgrp['time'].attrs['unit']

           
            destgrp['data'].attrs['unit'] = 'G'
            destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            
            
            del(bx,by,bz)
            del(xcal,ycal,zcal,xpol,ypol,zpol)
            if motion == 'fixed_rotation':
                del(x,y,z,rx,ry,rz,roll, pitch, yaw)
                
            if verbose:
                print("End of BDOT routine!")





if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_PL11B.hdf5") )
    tdiode_hdf = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode.hdf5") )
    dest = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_save_full.hdf5") )
    
    print('reading')
    util.mem()
    tstart = util.timeTest()
    full_filepath = bdot_raw_to_full(src, dest, tdiode_hdf=tdiode_hdf, grid=True, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
