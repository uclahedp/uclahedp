#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
bdot.py: BDOT analysis package
--> bdotRawToFull(src, dest, tdiode_hdf=None, grid=False, verbose=False)
    Takes in a source HDF5 file, integrates and calibrates signal based on
    metadata attributes of source HDF5. Optionally corrects for time offesets
    between shots using a timing diode source. Optionally outputs position 
    gridded data.
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




def bdotRawToFull(src, dest, tdiode_hdf=None, grid=False, verbose=False):
    """ Integrates bdot data, calibrates output using information about the probe.
        Corrects for probe angle based on which drive is being used.

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing bdot data
            
        dest: hdfPath object
            Path string to location processed bdot data should be written out

        tdiode_hdf:  hdfPath object
            Path to a raw hdf5 file containing tdiode data. T
            
        grid: Boolean
            If grid is true, output will be written in cartesian grid array
            format, eg. [nti, nx, ny, nz, nreps, nchan]. Otherwise, output will
            be in [nshots, nti, nchan] format


    Returns
    -------
       True (if executes to the end)
    """ 
    # ******
    # Load data from the raw HDF file
    # ******
    with h5py.File(src.file, 'r') as sf:
        
        #Get the datagroup
        srcgrp = sf[src.group]
        
        #Create dictionary of attributes
        attrs = hdftools.readAttrs(srcgrp)
        
        #Check for keys always required by this function
        req_keys = ['brd1','brd2','brd3', 'chan1','chan2', 'chan3', 
                    'xarea', 'yarea', 'zarea',
                    'xatten', 'yatten', 'zatten', 'gain',
                    'xpol', 'ypol', 'zpol', 'roll', 
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z',
                    'dt', 'nturns']
       
        
        if  'pos' in srcgrp:
            pos = srcgrp['pos'][:] #Read the entire array in
            req_keys = req_keys + ['motion_format']
            #If pos array exists, there are keywords required for that too.
            if srcgrp.attrs['motion_format'] == 'fixed_rotation':
                req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
            else:
                raise ValueError("Motion format unrecognized: " + str(srcgrp['pos'].attrs['motion_format']) )
            motion_format = srcgrp.attrs['motion_format']
        else:
            #If no position information is given, a single explicit position
            #is required. 
            req_keys = req_keys + ['xpos', 'ypos', 'zpos']
            grid = False #Can't grid data if there's no pos array!
            motion_format = None
            
            
        #Process the required keys, throwing an error if any cannot be found
        csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        

        #Extract the shape of the source data
        nshots, nti, nchan = srcgrp['data'].shape
        
        
        #If keywords specify gridded output, run postools.grid
        #shotgridind is a list of grid and repetition indices for each shot num
        #that tells the program where to store each data trace in the array
        if grid:
            shotgridind, xaxes, yaxes, zaxes = postools.gridShotIndList(pos, precision=.1)
            nx,ny,nz = len(xaxes), len(yaxes), len(zaxes)
            nreps =int( np.floor(nshots/(nx*ny*nz)))
            
        #If tdiode_hdf is set, run the tdiode processing routine
        #TODO: Maybe have a separate keyword to pass an already-processed
        #tdiode file here? Otherwise, this runs for each probe separaetly.
        if tdiode_hdf is not None:
            if verbose:
                print("Loading tdiode array from file.")
            #Get an array of all the t0 indices
            t0indarr = tdiode.calcT0ind(tdiode_hdf)
            #Get an array of all the good shots and badshots (indices)
            badshots, goodshots = tdiode.findBadShots(tdiode_hdf)
            #Replace any bad shots with a standin avg value
            #(Later should overwrite bad shots with good neighboring shots)
            t0indarr[badshots] = int(np.median(t0indarr[goodshots]))
            #We will remove up to max_t0shift indices from each array such that
            #the t0 indices all line up.
            min_t0ind = np.min(t0indarr[goodshots])
            max_t0shift = np.max(t0indarr[goodshots]) - min_t0ind
            #Compute new nti
            nti = nti - max_t0shift 
            
            
        if verbose:
            print("Opening destination HDF file")

        #Open the destination file
        #This exists WITHIN the open statement for the source file, so the
        #source file is open at the same time.
        with h5py.File(dest.file, 'a') as df:
            
            #Throw an error if this group already exists
            if dest.group is not '/' and dest.group in df.keys():
                raise hdftools.hdfGroupExists(dest)
            
            destgrp = df.require_group(dest.group)

            
            #Copy over attributes
            hdftools.copyAttrs(srcgrp, destgrp)

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

            
            
            atten = np.array([attrs['xatten'][0],attrs['yatten'][0],attrs['zatten'][0]])
            #Atten is assumed to be in dB. We could use the units to check this,
            #but it always is in dB and it's just a stupid source for errors.
            #Print this warning message if the units are different, just in case
            if attrs['xatten'][1].decode("utf-8") != 'dB':
                print("WARNING: ATTEN UNITS DO NOT MATCH dB")
                print(attrs['xatten'][1].decode("utf-8"))
                print("CONVERTING ANYWAY: CHECK YOUR UNITS!")

            #Convert atten to dB (if units set to dB)
            atten = np.power([10,10,10], atten/20.0) # Convert from decibels   
           
            # dt -> s
            dt = ( attrs['dt'][0]*u.Unit(attrs['dt'][1])).to(u.s).value

            # area : mm^2 -> m^2
            xarea = (attrs['xarea'][0]*u.Unit(attrs['xarea'][1])).to(u.m ** 2).value
            yarea = (attrs['yarea'][0]*u.Unit(attrs['yarea'][1])).to(u.m ** 2).value
            zarea = (attrs['zarea'][0]*u.Unit(attrs['zarea'][1])).to(u.m ** 2).value

            gain = attrs['gain'][0]
            nturns = attrs['nturns'][0]

            xcal = 1.0e4*dt*atten[0]/gain/(nturns*xarea)
            ycal = 1.0e4*dt*atten[1]/gain/(nturns*yarea)
            zcal = 1.0e4*dt*atten[2]/gain/(nturns*zarea)
            
            xpol = attrs['xpol'][0]
            ypol = attrs['ypol'][0]
            zpol = attrs['zpol'][0]
            
            
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
                    
                #Time remaining printout
                if verbose and i % nstepsreport == 0 and i > 0:
                    tremain = np.mean(tperstep)*(nshots - i)
                    print(str(i) + '/' + str(nshots) + ' complete, ' + 
                          util.timeFormat(tremain) + ' remaining' )


                #If a tdiode hdf was supplied, calculate the index correction
                #here
                if tdiode_hdf is not None:
                    ta = t0indarr[i] - min_t0ind
                    tb = ta + nti
                else:
                    #By default, read in the entire dataset
                    ta = 0
                    tb = -1
                
                #Read in the data from the source file
                bx = srcgrp['data'][i,ta:tb, 0]
                by = srcgrp['data'][i,ta:tb, 1]
                bz = srcgrp['data'][i,ta:tb, 2]
                
                #Apply the calibration factors
                bx = np.cumsum( detrend(bx) )*xcal*xpol
                by = np.cumsum( detrend(by) )*ycal*ypol
                bz = np.cumsum( detrend(bz) )*zcal*zpol
                
                #If a motion_format is set, apply the appropriate probe angle correction
                if motion_format == 'fixed_rotation':
                    #x,y,z is the probe's current position
                    x,y,z = srcgrp['pos'][i, :]
                    #rx, ry, rz is the location of the probe rotation point
                    #i.e. the center of the ball valve.
                    rx, ry, rz = attrs['rot_center_x'][0],attrs['rot_center_y'][0],attrs['rot_center_z'][0]
                    #x-rx, y-ry, z-rz is a vector pointing along the probe
                    #shaft towards the probe tip
                    #pitch is the angle of the probe shaft to the xz plane
                    pitch = np.arctan( (y-ry) / (x-rx) )
                    #yaw is the angle of the probe shaft to the xy plane
                    yaw = np.arctan( (z-rz) / (x-rx) ) 
                    
                    #If the probe is coming from the +X direction, its calibrated Z axis is already off by 180 degrees.
                    #This corrects for that, so the probes can be oriented wth Y+ up on both sides of LAPD
                    if (x-rx < 0.0):
                        yaw = yaw + np.pi
                    
                    #Roll is rotation of the probe about its axis, with
                    #y+ oriented up as roll=0
                    #This should be zero, unless a probe was later discovered
                    #to be incorrectly calibrated, so that the +Y mark was
                    #wrong
                    roll, unit = attrs['roll']
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
                    #Get location to write this datapoint from the shotgridind
                    xi = shotgridind[i, 0]
                    yi = shotgridind[i, 1]
                    zi = shotgridind[i, 2]
                    repi = shotgridind[i, 3]
                    #Write data
                    destgrp['data'][:, xi, yi, zi, repi, 0] = bx
                    destgrp['data'][:, xi, yi, zi, repi, 1] = by
                    destgrp['data'][:, xi, yi, zi, repi, 2] = bz
                else:
                    #Write data
                    destgrp['data'][i,:, 0] = bx
                    destgrp['data'][i,:, 1] = by 
                    destgrp['data'][i,:, 2] = bz                      
            
            
            
            
            if verbose:
                print("Writing axes to destination file")
            
            
            #Write the axes as required by the format of the data written
            if motion_format is not None:
                #Add the other axes and things we'd like in this file
                destgrp.require_dataset('pos', (nshots, 3), np.float32, chunks=True)[:] = srcgrp['pos'][:]
                for k in srcgrp['pos'].attrs.keys():
                    destgrp['pos'].attrs[k] = srcgrp['pos'].attrs[k]

            if grid:
                dimlabels = ['time', 'xaxes', 'yaxes', 'zaxes', 'reps', 'chan']
                
                destgrp.require_dataset('xaxes', (nx,), np.float32, chunks=True)[:] = xaxes
                destgrp['xaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']
                
                destgrp.require_dataset('yaxes', (ny,), np.float32, chunks=True)[:] = yaxes
                destgrp['yaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']
                
                destgrp.require_dataset('zaxes', (nz,), np.float32, chunks=True)[:] = zaxes
                destgrp['zaxes'].attrs['unit'] = srcgrp['pos'].attrs['unit']

            else:
                dimlabels = ['shots', 'time', 'chan']
                destgrp.require_dataset('shots', (nshots,), np.int32, chunks=True)[:] = srcgrp['shots'][:]
                destgrp['shots'].attrs['unit'] = srcgrp['shots'].attrs['unit']
                
            
            destgrp.require_dataset('chan', (nchan,), np.int32, chunks=True)[:] = srcgrp['chan'][:]
            destgrp['chan'].attrs['unit'] = srcgrp['chan'].attrs['unit']
            
            destgrp.require_dataset('time', (nti,), np.float32, chunks=True)
            destgrp['time'][:] = t
            destgrp['time'].attrs['unit'] = srcgrp['time'].attrs['unit']

           
            destgrp['data'].attrs['unit'] = 'G'
            destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            
            
            del(bx,by,bz)
            del(xcal,ycal,zcal,xpol,ypol,zpol)

            if verbose:
                print("End of BDOT routine!")
                
            return True





if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_LAPD1.hdf5") )
    tdiode_hdf = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode.hdf5") )
    dest = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_save_full.hdf5") )
    
    print('reading')
    util.mem()
    tstart = util.timeTest()
    full_filepath = bdotRawToFull(src, dest, tdiode_hdf=tdiode_hdf, grid=True, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
