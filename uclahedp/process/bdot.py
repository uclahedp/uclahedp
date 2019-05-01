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
import numpy as np
import os
import h5py
from scipy.signal import detrend as detrend
import astropy.units as u

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util
from uclahedp.tools import pos as postools
from uclahedp.tools import math



def bdotRawToFull(src, dest, tdiode_hdf=None, grid=False, 
                  verbose=False, debug = False,
                  offset_range=(0,100), offset_rel_t0 = (False, False), 
                  grid_precision=0.1, strict_grid=False, strict_axes = False):
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
            
        offset_range: tuple
            Tuple of indices between which the average of the signal will be
            computed and subtracted from the entire signal to correct for
            offset. This should be a segment with just noise, ideally at the
            very beginning of the dataset. Longer is better. 
            Default is (0,100)
            
        offset_rel_t0: Tuple of booleans
            If either of these values is set to True, the coorresponding
            offset_range value will be taken to be relative to the t0 index
            for that each shot. For example, if t0=2000 for a shot, 
            offset_range=(10, -100), and offset_rel_t0 = (False, True), then
            the offset will be computed over the range (10, 1900)
            
            
        grid_precision: float
            This is the precision to which position values will be rounded
            before being fit onto the grid. Only applies to fuzzy axis and grid
            creation.
            
        strict_axes: boolean
            If true, attempt to calculate axes from saved grid parameters.
            Default is false, which attempts to calculate axes by looking at
            position values.
            
        strict_grid: boolean
            If true, strictly unravel data onto the axes, assuming the probe
            moved in order reps->X->Y->Z. This will NOT correctly handle
            points where the probe was not at the requested position. Default
            is false, which applys "fuzzy gridding", which tries to find the
            best grid position for each shot individually.


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
        req_keys = ['xarea', 'yarea', 'zarea',
                    'xatten', 'yatten', 'zatten', 'gain',
                    'xpol', 'ypol', 'zpol', 'roll', 
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z',
                    'dt', 'nturns']
       


        if  'pos' in srcgrp:
            pos = srcgrp['pos'][:] #Read the entire array in
            #If pos array exists, there are keywords required for that too.
            motion_format = attrs['motion_format'][0]
            if motion_format == 'fixed_pivot':
                req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
            elif motion_format == 'cartesian':
                pass
            else:
                raise ValueError("Motion format unrecognized: " + str(attrs['motion_format'][0]) )
            
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
        
        #If requested by keyword, apply gridding
        if grid:
            
            #If grid parameters exist in the attrs, generate an exact grid
            #If not, generate a guess
            if strict_axes:
                try:
                    print("Applying stric axis creation")
                    nx = attrs['nx'][0]
                    ny = attrs['ny'][0]
                    nz = attrs['nz'][0]
                    dx = attrs['dx'][0]
                    dy = attrs['dy'][0]
                    dz = attrs['dz'][0]
                    x0 = attrs['x0'][0]
                    y0 = attrs['y0'][0]
                    z0 = attrs['z0'][0]
                    xaxis, yaxis, zaxis = postools.makeAxes(nx,ny,nz,
                                                            dx,dy,dz,
                                                            x0,y0,z0)
                except KeyError:
                    print("Missing axis parameters: attempting fuzzy axis creation")
                    xaxis,yaxis,zaxis = postools.guessAxes(pos, precision=grid_precision)
            else:
                print("Applying fuzzy axis creation")
                xaxis,yaxis,zaxis = postools.guessAxes(pos, precision=grid_precision)
                
            #Calculate length of axes
            nx, ny, nz, nreps = postools.calcNpoints(pos, xaxis, yaxis, zaxis)
            
            #This line SHOULD be redundent, but sometimes it's necessary.
            #It is possible, when combining two motion lists, to get
            #some extra shots at the end of the datarun that don't have
            #positions associated. Recalculating this here ensures we only
            #take the ones relevant to this grid
            nshots = nx*ny*nz*nreps
            
            #If grid precision is zero, apply strict gridding
            #Otherwise, apply fuzzy gridding
            if strict_grid:
                print("Applying strict gridding")
                shotgridind = postools.strictGrid(nx,ny,nz,nreps)  
            else:
                print("Applying fuzzy gridding")
                shotgridind = postools.fuzzyGrid(pos, xaxis, yaxis, zaxis, 
                                                 precision=grid_precision)
            
            if debug:
                print("** Gridding parameters **")
                print('nshots: ' + str(nshots))
                print('nx, ny, nz, nreps: ' + str((nx, ny, nz, nreps)))
                print('dx, dy, dz: ' + str((dx, dy, dz)))
                print('x0, y0, z0: ' + str((x0, y0, z0)))

            
            
        #If tdiode_hdf is set, load the pre-processed tdiode data
        if tdiode_hdf is not None:
            if verbose:
                print("Loading tdiode array from file.")
            with h5py.File(tdiode_hdf.file, 'r') as sf:
                grp = sf[tdiode_hdf.group]
                t0indarr = grp['t0indarr'][:]
                goodshots = grp['goodshots'][:]
            #We will remove up to max_t0shift indices from each array such that
            #the t0 indices all line up.
            min_t0ind = np.min(t0indarr[goodshots])
            max_t0shift = np.max(t0indarr[goodshots]) - min_t0ind
            #Compute new nti
            print(max_t0shift)
            nti = nti - max_t0shift 
            
        if verbose:
            print("Opening destination HDF file")
        
        #Create the destination file directory if necessary
        hdftools.requireDirs(dest.file)

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
                destgrp.require_dataset('data', (nti, nx, ny, nz, nreps, nchan), np.float32, chunks=(np.min([nti, 20000]),1,1,1,1,1), compression='gzip')
            else:
                destgrp.require_dataset('data', (nshots, nti, nchan), np.float32, chunks=(1, np.min([nti, 20000]), 1), compression='gzip')
            
            #Load the time vector
            t = srcgrp['time']
            #If a timing diode is being applied, correct the time vector here.
            if tdiode_hdf is not None:
                t = t[0:nti] - t[min_t0ind]
    
            #Asssemble the array of calibration factors from the attrs dict
            cal = calibrationFactor(attrs)
            
            

            #Initialize time-remaining printout
            tr = util.timeRemaining(nshots)
            
            if verbose:
                print("Beginning processing data shot-by-shot.")
            
            #Chunking data processing loop limits memory usage
            for i in range(nshots):
                
                #Update time remaining
                if verbose:
                        tr.updateTimeRemaining(i)

                #If a tdiode hdf was supplied, calculate the index correction
                #here
                if tdiode_hdf is not None:
                    #Calculate the starting and ending arrays for the data
                    ta = t0indarr[i] - min_t0ind
                    tb = ta + nti

                    #Calculate the range over which to calculate the offset
                    #for each shot
                    #If offset_rel_t0 is set for either point, add the t0 array
                    if offset_rel_t0[0]:
                        offset_a = offset_range[0] + t0indarr[i]
                    else:
                        offset_a = offset_range[0]
                        
                    if offset_rel_t0[1]:
                        offset_b = offset_range[1] + t0indarr[i]
                    else:
                        offset_b = offset_range[1]
                    
                else:
                    #By default, read in the entire dataset
                    ta = None
                    tb = None
                    offset_a = offset_range[0]
                    offset_b = offset_range[1]
                    
                if debug:
                    print("Data range: [" + str(ta) + "," + str(tb) + "]")
                    print("Offset range: [" + str(offset_a) + "," + 
                                          str(offset_b) + "]")
                    
                    
                
                #Read in the data from the source file
                bx = srcgrp['data'][i,ta:tb, 0]
                by = srcgrp['data'][i,ta:tb, 1]
                bz = srcgrp['data'][i,ta:tb, 2]
                
                
                #Remove offset from each channel
                bx = bx - np.mean(bx[offset_a:offset_b])
                by = by - np.mean(by[offset_a:offset_b])
                bz = bz - np.mean(bz[offset_a:offset_b])
                
                #Apply the calibration factors
                bx = np.cumsum( bx )*cal[0]
                by = np.cumsum( by )*cal[1]
                bz = np.cumsum( bz )*cal[2]
                
                #If a motion_format is set, apply the appropriate probe angle correction
                if motion_format == 'cartesian':
                    #Don't need to make any correction
                    pass 
                elif motion_format == 'fixed_pivot':
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
                    
                    #If the probe is coming from the -X direction, its calibrated Z axis is already off by 180 degrees.
                    #This is because the probes are calibrated to match the East side of LAPD
                    if ((x-rx) > 0.0):
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
                    try:
                        destgrp['data'][:, xi, yi, zi, repi, 0] = bx
                        destgrp['data'][:, xi, yi, zi, repi, 1] = by
                        destgrp['data'][:, xi, yi, zi, repi, 2] = bz
                    except ValueError as e:
                        print("ERROR!")
                        print(destgrp['data'].shape)
                        print(bx.shape)
                        print([xi, yi, zi, repi])
                        raise(e)
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
                destgrp.require_dataset('pos', (nshots, 3), np.float32, chunks=True)[:] = srcgrp['pos'][0:nshots]
                for k in srcgrp['pos'].attrs.keys():
                    destgrp['pos'].attrs[k] = srcgrp['pos'].attrs[k]

            if grid:
                dimlabels = ['time', 'xaxis', 'yaxis', 'zaxis', 'reps', 'chan']
                
                destgrp.require_dataset('xaxis', (nx,), np.float32, chunks=True)[:] = xaxis
                destgrp['xaxis'].attrs['unit'] = srcgrp.attrs['motion_unit']
                
                destgrp.require_dataset('yaxis', (ny,), np.float32, chunks=True)[:] = yaxis
                destgrp['yaxis'].attrs['unit'] = srcgrp.attrs['motion_unit']
                
                destgrp.require_dataset('zaxis', (nz,), np.float32, chunks=True)[:] = zaxis
                destgrp['zaxis'].attrs['unit'] = srcgrp.attrs['motion_unit']
                
                destgrp.require_dataset('reps', (nreps,), np.int32, chunks=True)[:] = np.arange(nreps)
                destgrp['reps'].attrs['unit'] = ''

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

            if verbose:
                print("End of BDOT routine!")
                
            return True
        
        
def calibrationFactor(attrs):
    atten = np.array([attrs['xatten'][0],attrs['yatten'][0],attrs['zatten'][0]])
    #Atten is assumed to be in dB. We could use the units to check this,
    #but it always is in dB and it's just a stupid source for errors.
    #Print this warning message if the units are different, just in case
    if attrs['xatten'][1] != 'dB':
        print("WARNING: ATTEN UNITS DO NOT MATCH dB")
        print(attrs['xatten'][1])
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
    
    xpol = attrs['xpol'][0]
    ypol = attrs['ypol'][0]
    zpol = attrs['zpol'][0]

    xcal = 1.0e4*dt*atten[0]/gain/(nturns*xarea)*xpol
    ycal = 1.0e4*dt*atten[1]/gain/(nturns*yarea)*ypol
    zcal = 1.0e4*dt*atten[2]/gain/(nturns*zarea)*zpol
    
    cal = [xcal, ycal, zcal]
    
    return cal




def fullToCurrent(src, dest, verbose=False):
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        try:
            dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
            shape =  srcgrp['data'].shape
        except KeyError: 
            raise KeyError("bdot.fullToCurrent requires the data array to have an attribute 'dimensions' and 'shape'")
            
        #We will duplicate the chunking on the new array
        chunks = srcgrp['data'].chunks
        

        try:
            xax = dimlabels.index("xaxis") 
            yax = dimlabels.index("yaxis") 
            zax = dimlabels.index("zaxis") 
            
            xaxis = srcgrp['xaxis']
            yaxis = srcgrp['yaxis']
            zaxis = srcgrp['zaxis']
            
            nti = shape[ dimlabels.index("time")  ]
            nx = shape[xax]
            ny = shape[yax]
            nz = shape[zax]
            
        except KeyError:
            raise KeyError("bdot.fullToCurrent requires dimensions 'time', 'xaxis', 'yaxis', 'zaxis'")
            
            
        if nti > 10000:
            print("WARNING: NTI IS LARGE! CURRENT CALCULATION WILL TAKE A VERY LONG TIME!")
            print("If you have better things to do with your CPU hours, try thinning the data first.")
        
        #Create the destination file directory if necessary
        hdftools.requireDirs(dest.file)
        
        with h5py.File(dest.file, 'w') as df:
            destgrp = df[dest.group]
            
            destgrp.require_dataset('data', shape, np.float32, chunks=chunks, compression='gzip')
            destgrp['data'].attrs['unit'] = 'A/cm^2'
            destgrp['data'].attrs['dimensions'] = hdftools.strListToArr(dimlabels)
            
            #Copy the axes over
            for ax in dimlabels:
                srcgrp.copy(ax, destgrp)
                
                
            chunksize = 100
            nchunks = int(np.ceil(nti/chunksize))
            
            #Initialize time-remaining printout
            tr = util.timeRemaining(nchunks)
            
            for i in range(nchunks):
                #Update time remaining
                if verbose:
                        tr.updateTimeRemaining(i)

                a = i*chunksize
                if i == nchunks-1:
                    b = None
                else:
                    b = (i+1)*chunksize
                
                #Constant is (c/4pi) * (conversion CGS -> A/m^2)*(conversion A/m^2 -> A/cm^2)
                #(2.99e10/4pi)*(3.0e-5)*(1e-4)
                #3e-5 is from the NRL formulary
                destgrp['data'][a:b, ...] = (7.138)*math.curl(srcgrp['data'][a:b, ...], 
                    xax, yax, zax, xaxis, yaxis, zaxis)
                
        return dest

    




if __name__ == "__main__":
    #raw = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_raw.hdf5") )
    #tdiode_hdf = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_tdiode_raw.hdf5") )
    #full = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_full.hdf5") )
    #current = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_current.hdf5") )
    
    exp = 'LAPD_Jan2019'
    probe = 'LAPD_C6'
    run = 30
    
    src = hdftools.hdfPath( os.path.join("F:", exp, "RAW", 'run' + str(run) + '_' + probe + '_raw.hdf5'))
    tdiode_hdf = hdftools.hdfPath(os.path.join("F:", exp, "FULL", 'run' + str(run) + '_' + 'tdiode' + '_full.hdf5'))
    dest = hdftools.hdfPath(os.path.join("F:", exp, "FULL", 'run' + str(run) + '_' + probe + '_full.hdf5'))
    
    src = hdftools.hdfPath( '/Volumes/PVH_DATA/' + exp + '/RAW/run' + str(run) + '_' + probe + '_raw.hdf5')
    tdiode_hdf = hdftools.hdfPath('/Volumes/PVH_DATA/' + exp + '/FULL/' + 'run' + str(run) + '_' + 'tdiode' + '_full.hdf5')
    dest = hdftools.hdfPath('/Volumes/PVH_DATA/'+ exp + '/FULL/' + 'run' + str(run) + '_' + probe + '_full.hdf5')
    

    #Delete the output file if it already exists
    try:
        os.remove(dest.file)
    except FileNotFoundError:
        pass
    
    print('reading')
    util.mem()
    tstart = util.timeTest()
    print(src.file)
    full_filepath = bdotRawToFull(src, dest, tdiode_hdf=tdiode_hdf, grid=True, verbose=True, debug=False, 
                                  offset_range = (0, -100), offset_rel_t0 = (False, True), 
                                  strict_axes = True, strict_grid = False, grid_precision=0.1)
    #cur_filepath = fullToCurrent(dest, current, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
