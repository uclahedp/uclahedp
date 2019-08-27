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

from scipy.optimize import curve_fit as curve_fit
import csv
import matplotlib.pyplot as plt

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util
from uclahedp.tools import pos as postools
from uclahedp.tools import math



def bdotRawToFull(src, dest, 
                  tdiode_hdf=None, grid=False, integrate=True, 
                  calibrate =True, highfreq_calibrate=True,
                  angle_correction = True, remove_offset = True,
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
            Path to a raw hdf5 file containing tdiode data. If no HDF file is
            provided, no timing correction will be applied.
            
        grid: Boolean
            If grid is true, output will be written in cartesian grid array
            format, eg. [nti, nx, ny, nz, nreps, nchan]. Otherwise, output will
            be in [nshots, nti, nchan] format
            
            
        integrate: Boolean
             If True, integrate the bdot data (usually you want to do this).
             Default is True
             
        calibrate: Boolean
             If True, calculate and apply ANY calibration factors 
             to the data. Default is True.
             
        highfreq_calibrate: Boolean
             If True, calculate and apply the high frequency calibration 
             factors to the data. Default is True. If the 'tau' variables are
             not specified in the probe metadata, the HF calibration won't be
             applied regardless of this keyword.
             
        angle_correction: Boolean
             If True, apply any angular correction between axes that is
             required based on the motion_format keyword in the metadata. If 
             false, no correction is applied regardless of the metadata.
             Default is True.
             
       remove_offset: Boolean
            If True, remove an offset from the data based on the offset_range
            specified in those keywords. If False, data will remain as-is. 
            Default is True.
            
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
            if motion_format == 'fixed_pivot' and angle_correction:
                req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
            elif motion_format == 'cartesian' and angle_correction:
                pass
            elif not angle_correction:
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
           shotgridind, xaxis, yaxis, zaxis, nx, ny, nz, nreps, nshots = postools.grid(
                     pos, attrs, strict_axes=strict_axes, 
                     strict_grid=strict_grid, grid_precision=grid_precision, 
                     invert=False)
          
            
        

            
            
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
    
          
            if calibrate:
                 # dt -> s
                 dt = ( attrs['dt'][0]*u.Unit(attrs['dt'][1])).to(u.s).value
                
                 #First calculate the low frequency calibration factors
                 calAx, calAy, calAz = calibrationFactorsLF(attrs)
                 
                 #If HF calibration factors are provided, calculate those
                 #calibraton constants too
                 if 'xtau' in attrs.keys() and highfreq_calibrate:
                     calBx, calBy, calBz = calibrationFactorsHF(attrs)
                 else:
                      calBx, calBy, calBz = None,None,None
                      
                 print(calBx)
            
            

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
                if tdiode_hdf is not None and remove_offset:
                    #Calculate the starting and ending arrays for the data
                    ta = t0indarr[i] - min_t0ind
                    tb = ta + nti

                    #Calculate the range over which to calculate the offset
                    #for each shot
                    #If offset_rel_t0 is set for either point, add the t0 array
                    if offset_rel_t0[0]:
                        offset_a = offset_range[0] + t0indarr[i] - ta
                    else:
                        offset_a = offset_range[0]
                        
                    if offset_rel_t0[1]:
                        offset_b = offset_range[1] + t0indarr[i] - ta
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
                dbx = srcgrp['data'][i,ta:tb, 0]
                dby = srcgrp['data'][i,ta:tb, 1]
                dbz = srcgrp['data'][i,ta:tb, 2]
                
                
                if remove_offset:
                     #Remove offset from each channel
                     dbx = dbx - np.mean(dbx[offset_a:offset_b])
                     dby = dby - np.mean(dby[offset_a:offset_b])
                     dbz = dbz - np.mean(dbz[offset_a:offset_b])
                     
                                #Integration comes after calibration?
                if integrate:
                     #Intgrate
                     bx = np.cumsum(dbx)*dt
                     by = np.cumsum(dby)*dt
                     bz = np.cumsum(dbz)*dt
                else:
                    bx,by,bz = dbx, dby, dbz
                
                
                if calibrate:
                     #Apply the high-frequency calibration if one was
                     #provided
                     if calBx is not None and highfreq_calibrate:
                          bx = bx + calBx*dbx
                          by = by + calBy*dby
                          bz = bz + calBz*dbz

                     #Apply the low-frequency calibration factors
                     #Probe pol dir is included in these
                     bx = bx*calAx
                     by = by*calAy
                     bz = bz*calAz
                
                
                #If a motion_format is set, apply the appropriate probe angle correction
                if motion_format == 'cartesian' and angle_correction:
                    #Don't need to make any correction
                    pass 
                elif motion_format == 'fixed_pivot' and angle_correction:
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
                destgrp['xaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
                destgrp.require_dataset('yaxis', (ny,), np.float32, chunks=True)[:] = yaxis
                destgrp['yaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
                destgrp.require_dataset('zaxis', (nz,), np.float32, chunks=True)[:] = zaxis
                destgrp['zaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
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

           
            if calibrate:
                 destgrp['data'].attrs['unit'] = 'G'
            else:
                 destgrp['data'].attrs['unit'] = 'V'
                 
            destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            
            
            del(bx,by,bz)

            if verbose:
                print("End of BDOT routine!")
                
            return True
        

def calibrationFactorsLF(attrs):
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

    # area : mm^2 -> m^2
    xarea = (attrs['xarea'][0]*u.Unit(attrs['xarea'][1])).to(u.m ** 2).value
    yarea = (attrs['yarea'][0]*u.Unit(attrs['yarea'][1])).to(u.m ** 2).value
    zarea = (attrs['zarea'][0]*u.Unit(attrs['zarea'][1])).to(u.m ** 2).value

    gain = attrs['gain'][0]
    nturns = attrs['nturns'][0]
    
    xpol = attrs['xpol'][0]
    ypol = attrs['ypol'][0]
    zpol = attrs['zpol'][0]

    xcal = 1.0e4*atten[0]/gain/(nturns*xarea)*xpol
    ycal = 1.0e4*atten[1]/gain/(nturns*yarea)*ypol
    zcal = 1.0e4*atten[2]/gain/(nturns*zarea)*zpol
    
    return xcal, ycal, zcal


def calibrationFactorsHF(attrs):
    # area : convert to seconds
    xtau = (attrs['xtau'][0]*u.Unit(attrs['xtau'][1])).to(u.s).value
    ytau = (attrs['ytau'][0]*u.Unit(attrs['ytau'][1])).to(u.s).value
    ztau = (attrs['ztau'][0]*u.Unit(attrs['ztau'][1])).to(u.s).value
    
    calBx = xtau
    calBy = ytau
    calBz = ztau

    return calBx, calBy, calBz



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
            tr = util.timeRemaining(nchunks, reportevery=10)
            
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

#These two functions are used in the calibrateProbe routine
def lfProbeArea(freq, mag, nturns, hturns, gain,  Rp, r):
     mu0 = 4*np.pi*1e-7
     return mag*Rp*r/(gain*hturns*np.power(4/5, 1.5)*mu0*nturns)/(2*np.pi*freq)


def lfCoil(freq, nturns, hturns, gain, area, Rp, r):
     
     mu0 = 4*np.pi*1e-7
     w = 2*np.pi*freq
     
     coeff = np.power(4/5, 1.5)*hturns*mu0*area*nturns*gain/(r*Rp)
     
     #This function just returns the imaginary part, since that's what is used
     #real included here for completeness...
     #real = coeff*np.power(w,2)*tau*tdelay
     im = coeff*w
     return im

def hfCoil(freq, nturns, hturns, gain, area, Rp, r, tau, tdelay):
     mu0 = 4*np.pi*1e-7
     
     w = 2*np.pi*freq

     coeff = np.power(4/5, 1.5)*hturns*mu0*area*nturns*gain/(r*Rp)
     x = w/(1 + np.power(w*tau,2))
     A = w*tau*np.cos(w*tdelay) - np.sin(w*tdelay)
     B = w*tau*np.sin(w*tdelay) + np.cos(w*tdelay)
     #return area*nturns*gain*(16/np.power(5,1.5))*mu0/(r*Rp)*(tau-tdelay)*np.power(freq,2)
     real = coeff*x*A
     im = coeff*x*B

     return np.concatenate((real,im))
    
def calibrateProbe(file, nturns, gain, hturns=32, Rp=10, r=0.055, area_freq_range = [1e2,1e6]):
     """
     csvfile -> Bdot calibration csv file with the following columns...
     
     "freq" -> Frequencies in Hz
     
     For s in [x,y,z] EITHER
     "smag" and "sphase" with magnitude in dB and phase in degrees
     OR
     "sreal" and "sim" with real and imaginary parts in dB
     
     
     Required Keywords (Likely to change)
     nturns -> number of bdot coil turns
     gain -> gain of differential amplifier used for calibration
     
     Other Keywords (less likely to change bc they are part of the LAPD's bdot
     testing setup)
     
     Rp -> Resistance of the resistor used for measuring the coil current.
     r -> Radius of the Helmholtz coil in meters
     
     
     """

     
     #Determine the type of file being supplied
     ext = os.path.splitext(file)[1].lower()


     if ext == '.dat':
         print("Binary files not currently supported here: ask Pat for his" +
               " converted script.")

    
     elif ext == '.csv':
         #Figure out the number of lines in the file, assuming one header row
         #There really doesn't seem to be a better way of doing this?
         with open(file) as csvfile:
              reader = csv.DictReader(csvfile)
              nlines = sum(1 for row in reader) - 1

         freq = np.zeros(nlines)
         #Signal contains mag, phase, real, imaginary in that order for each 
         # of the three channels
         sig = np.zeros([nlines, 4, 3])
         lf_fit = np.zeros([nlines,3])
         hf_fit_real = np.zeros([nlines,3])
         hf_fit_im = np.zeros([nlines,3])

         #Read the file as a dictionary
         with open(file) as csvfile:
              reader = csv.DictReader(csvfile)
              keys = reader.fieldnames
              
              #Determine whether the file contains magnitude/phase data
              #or real/imaginary data
              #The network analyzer stores both...
              if 'xmag' in keys:
                   mag_phase = True
              elif 'xreal' in keys:
                   mag_phase = False
              else:
                   raise(KeyError("No mag/phase or real/im keys found!"))
                   
              #Adjust index because of header
              header = next(reader)
              for i, row in enumerate(reader):
                   #This is assumed to be in Hz
                   freq[i] = float(row['freq'])
                   
                   if mag_phase:
                        #Magnitudes are all assumed to be in dB
                        sig[i, 0, 0] = pow(10.0, float(row['xmag'])/20.0)
                        sig[i, 0, 1] = pow(10.0, float(row['ymag'])/20.0)
                        sig[i, 0, 2] = pow(10.0, float(row['zmag'])/20.0)
                   
                        #Phase is assumed to be in degrees
                        #Why ths network analyzer does this...who knows
                        sig[i, 1, 0] = np.radians(  float(row['xphase'])  )
                        sig[i, 1, 1] = np.radians(  float(row['xphase'])  )
                        sig[i, 1, 2] = np.radians(  float(row['xphase'])  )
         
                   else:
                        #TODO: test this part? (I don't have a file like this handy)
                        #Real and imaginary parts are all assumed to be in dB
                        sig[i, 2, 0] = pow(10.0, float(row['xreal'])/20.0)
                        sig[i, 2, 1] = pow(10.0, float(row['yreal'])/20.0)
                        sig[i, 2, 2] = pow(10.0, float(row['zreal'])/20.0)
                        sig[i, 3, 0] = pow(10.0, float(row['xim'])/20.0)
                        sig[i, 3, 1] = pow(10.0, float(row['yim'])/20.0)
                        sig[i, 3, 2] = pow(10.0, float(row['zim'])/20.0)
          
              
        
     else:
         print("File extension not supported:" + ext)
         

     #Whichever data wasn't read in, calculate it          
     if mag_phase:
        #real = mag*cos(phase)
        sig[:,2,:] = sig[:,0,:]*np.cos(sig[:,1,:])
        #imaginary = mag*sin(phase)
        sig[:,3,:] = sig[:,0,:]*np.sin(sig[:,1,:])
     else:
        #mag = sqrt(real^2 + im^2)
        sig[:,0,:] = np.sqrt(np.power(sig[:,2,:],2)  + np.power(sig[:,3,:],2))
        #phase = arctan(im/real)
        sig[:,1,:] = np.arctan(sig[:,3,:]/sig[:,2,:])
       

     area = np.zeros(3)
     tau = np.zeros(3)
     tdelay = np.zeros(3)
      
      
     a = np.argmin(np.abs(freq - area_freq_range[0]))
     b = np.argmin(np.abs(freq - area_freq_range[1] ))
      
     for i in range(3):
          #Calculate the area of the coil data from just the specified
          #frequency range
          #area[i] = np.median( lfProbeArea(freq[a:b], sig[a:b,0,i], 
          #                        nturns, hturns,gain, Rp, r))
           
          fcn = lambda freq, area: lfCoil(freq,
                                          nturns, hturns, gain, area, 
                                          Rp, r)
           
          #Fit over the range a:b to find the area
          popt, pcov = curve_fit(fcn, freq[a:b], sig[a:b,3,i], 
                                                p0=[ 1])
          area[i] = popt[0]
          
          lf_fit[:,i] = fcn(freq,  area[i])
           
           
          #Now fit the full signal (with the area fixed) for the impedence 
          #time constant tau
          fcn = lambda freq, tau, tdelay: hfCoil(freq, 
                                                 nturns, hturns, gain, 
                                                 area[i], 
                                                 Rp, r, tau, tdelay) 
           
          concat_data = np.concatenate((sig[:,2,i],sig[:,3,i]))
          
          popt, pcov = curve_fit(fcn, freq, concat_data, 
                                                 p0=[ 1e-7,-1e-7])
          tau[i] = popt[0]
          tdelay[i] = popt[1]
          
          hf_fit_real[:,i] = fcn(freq,  tau[i], tdelay[i])[0:nlines]
          hf_fit_im[:,i] = fcn(freq,  tau[i], tdelay[i])[nlines:2*nlines]
 

     print("**** Bdot Calibration Report *****")
     print("File:" + str(file))
     for i in range(3):
          axes = ['x', 'y', 'z']
          print("************")
          print(axes[i] + 'area: ' + str(np.round(area[i]*1e6, decimals=3)) + ' mm^2')
          print(axes[i] + 'tau: ' + str(np.round(tau[i]*1e9, decimals=3)) + ' ns')
          print(axes[i] + 'tdelay: ' + str(np.round(tdelay[i]*1e9, decimals=3)) + ' ns')
           
           
     fig, ax = plt.subplots( nrows=3, ncols=3, figsize = [8,8])
     fig.subplots_adjust( hspace=.35, wspace=.35)
     fontsize = 12
     for i in range(3):
          lf_xrange = (area_freq_range[0]*1e-3,area_freq_range[1]*1e-3)
          hf_xrange = (np.min(freq)*1e-3,np.max(freq)*1e-3)
          #hf_xrange = (20, 5000)

          if i ==2:
               ax[i,0].set_xlabel('Frequency (kHz)', fontsize=fontsize)
               ax[i,1].set_xlabel('Frequency (kHz)', fontsize=fontsize)
               ax[i,2].set_xlabel('Frequency (kHz)', fontsize=fontsize)
          ax[i,0].set_title(axes[i] + ' Low Freq.')
          ax[i,0].plot(freq[a:b]*1e-3, sig[a:b:,3,i], linewidth=3)
          ax[i,0].plot(freq[a:b]*1e-3, lf_fit[a:b,i])
          ax[i,0].set_xlim(lf_xrange)
          ax[i,0].set_ylabel('Im( V$_{m}$ / V$_{o}$ )', fontsize=fontsize)
           
          ax[i,1].set_title(axes[i] + ' High Freq. Real')
          ax[i,1].plot(freq*1e-3, sig[:,2,i], linewidth=3)
          ax[i,1].plot(freq*1e-3, hf_fit_real[0:nlines,i])
          ax[i,1].set_xlim(hf_xrange)
          ax[i,1].set_xscale('log')
          ax[i,1].set_ylabel('Re( V$_{m}$ / V$_{o}$ )', fontsize=fontsize)
           
          ax[i,2].set_title(axes[i] + ' High Freq. Imaginary')
          ax[i,2].plot(freq*1e-3, sig[:,3,i], linewidth=3)
          ax[i,2].plot(freq*1e-3, hf_fit_im[0:nlines,i])
          ax[i,2].set_xlim(hf_xrange)
          ax[i,2].set_xscale('log')
          ax[i,2].set_ylabel('Im( V$_{m}$ / V$_{o}$ )', fontsize=fontsize)

               

if __name__ == "__main__":
     
     #csvfile = os.path.join("G:","LAPD_Mar2018","Bdot Calibration Data", "LAPD7.csv")
     #csvfile = os.path.join("/Volumes","PVH_DATA","LAPD_Mar2018","Bdot Calibration Data", "LAPD7.csv")
     #csvfile = os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","bdot_calibration", "LAPD_C2_BX.dat")
     #csvfile = os.path.join("G:","LAPD_Jul2019","bdot_calibration", "LAPD_C2_jeff.csv")
     
     #calibrateProbe(csvfile, 10, 100)
     
     
     src = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "run34_LAPD_C2_full.hdf5"))
     dest = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "run34_LAPD_C2_current.hdf5"))
     fullToCurrent(src, dest, verbose=False)
     
     
     """
    #raw = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_raw.hdf5") )
    #tdiode_hdf = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_tdiode_raw.hdf5") )
    #full = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_full.hdf5") )
    #current = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "run103_PL11B_current.hdf5") )
    
    exp = 'LAPD_Jan2019'
    probe = 'LAPD_C6'
    run = 25
    
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
    full_filepath = bdotRawToFull(src, dest, tdiode_hdf=tdiode_hdf, grid=True, verbose=True, debug=False, 
                                  offset_range = (0, -100), offset_rel_t0 = (False, True), 
                                  strict_axes = True, strict_grid = False, grid_precision=0.1)
    #cur_filepath = fullToCurrent(dest, current, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
    """
