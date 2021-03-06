# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 08:10:47 2019

@author: Peter
"""

import os, h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as curve_fit
from scipy.signal import find_peaks as find_peaks
import astropy.units as u

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util
from uclahedp.tools import pos as postools
from uclahedp.tools import math

#Lots of fitting warnings are thrown: this is used to catch and hide them
import warnings

def fixRange(array, bound):
     """
     
     """
     if bound is None:
          ind = np.shape(array)[0] - 1
     else:
          ind = np.argmin(np.abs(array - bound))
     return ind


def expFcn(V, A, B, kT):
     return A*np.exp((V)/kT) - B
     

def normalize(f):
     offset = np.median(f[0:100])
     f = f - offset
     factor = 1/np.max(f)
     f = f*factor
     return f, offset, factor


def choose_sweep(time, voltage, desired_time, plots=False):
     peaktimes, a, b = find_sweeps(time, voltage, plots=plots)
     ind = np.argmin(np.abs(peaktimes - desired_time))
     return peaktimes[ind], int(a[ind]), int(b[ind])

def find_sweeps(time, voltage, plots=False):
     #Normalize
     voltage, voffset, vfactor = normalize(voltage)
   
     if plots:
          fig_full, ax_full = plt.subplots(figsize = [8,4])
          fig_full.suptitle("Entire Ramp Signal")
          ax_full.set_xlabel("Time (ms)")
          ax_full.set_ylabel("Voltage Sweep (V)")
          ax_full.plot(time*1e3, voltage, 'k')
         

          
     #Define a mainumum height required to be a "peak"
     req_height = 0.9
     #Set a minimum spacing between "peaks" so only one index is returned for
     #each ramp
     req_separation = 0.01*len(time)
     #Find the peaks
     peaks, props = find_peaks(voltage, height=req_height, 
                               distance = req_separation)
     
     npeaks = peaks.size
 
     #Mark the peaks on the plot
     if plots:
         for p in peaks:
             ax_full.axvline(x=time[p]*1e3)
     
     peaktime = np.zeros(npeaks)
     start = np.zeros(npeaks)
     end = np.zeros(npeaks)
     
     for i,peak in enumerate(peaks):
          b = peaks[i]

          slope = np.mean(np.gradient(voltage[b-100:b]))
          a = b - int(1/slope)
          
          #print(str(a) + ':' + str(b))
          
  
          if plots:
               ax_full.plot(time[a:b]*1e3, voltage[a:b], color='blue', linewidth=1)
          
          peaktime[i] = time[int(np.mean([a,b]))]
          start[i] = int(a)
          end[i] = int(b)
          
     return peaktime, start, end
          
          
def vsweep_fit(voltage, current, esat_range=None, exp_range=None, plots=False, 
               verbose=False, return_fits=False, fail_val=0, area=1 ):
     """
     Fits several regions of the IV curve to determine plasma parameters, then
     calculates several others and returns all of them. 
     
     The fitting process is often very messy and generates many errors if the 
     data is not good. The algorithm deals with this by attempting to fit the
     data assuming it is good, and handling any errors thrown by replacing the
     data with a placeholder value fail_val.
     
     area: float (cm^2)
     Area of the probe tip
     
     """
     
     with warnings.catch_warnings():
         #Ignore fitting errors: there will be tons on bad data
         warnings.filterwarnings('ignore')
     
         #The entire fitting program is wrapped in a try/catch statement
         #If any error is thrown, all results will be reported as a dummy value "fail_val"
         #Some checks are made manually and will throw a ValueError if they fail, resulting
         #in the same dummy values being returned
         try:
    
             #This is a convention. I guess I'll follow it.
             current = -current
             
             if plots:
                  fig_input, ax_input = plt.subplots(ncols=2, figsize = [8,3])
                  fig_input.suptitle("Input Signal Summary")
                  fig_input.subplots_adjust( wspace=.4)
                  ax_input[0].set_title("Ramp")
                  ax_input[0].set_xlabel("Time indices")
                  ax_input[0].set_ylabel("Sweep Voltage (V)")
                  ax_input[0].plot(voltage, 'k')
             
                  ax_input[1].set_title("Langmuir Curve")
                  ax_input[1].set_xlabel("Sweep Voltage (V)")
                  ax_input[1].set_ylabel("Electron Current (A)")
                  ax_input[1].plot(voltage, current, 'k')
                  
             if esat_range is None or exp_range is None:
                  norm_i, i_offset, i_factor = normalize(current)
        
                  esat_ceil = 1.0
                  esat_floor = .9
                  
                  exp_ceil = 0.5
                  exp_floor = .05
                  
                  if esat_range is None:
                       esat_a = np.argmin(np.abs(esat_floor  - norm_i))
                       esat_b = np.argmin(np.abs(esat_ceil  - norm_i))
                       esat_range = [voltage[esat_a], voltage[esat_b]]
                  if exp_range is None:
                       exp_a = np.argmin(np.abs(exp_floor - norm_i))
                       exp_b = np.argmin(np.abs(exp_ceil  - norm_i))
                       exp_range = [voltage[exp_a], voltage[exp_b]]
                  
                  if plots:
                       fig_guess, ax_guess = plt.subplots(figsize = [8,3])
                       ax_guess.plot(voltage, norm_i, 'k')
                       ax_guess.axhspan(esat_floor, esat_ceil, alpha=0.25, color='green')
                       ax_guess.axhspan(exp_floor, exp_ceil, alpha=0.25, color='red')
        
             #Calculate indices based on the esat range
             esat_a = fixRange(voltage, esat_range[0])
             esat_b = fixRange(voltage, esat_range[1])
             
             #This indicates a bad curve that will not fit
             if esat_a > esat_b:
                 raise ValueError
             
             #Pull those parts of the curve
             esat_v = voltage[esat_a:esat_b]
             esat_i = current[esat_a:esat_b]
             #Fit and store fit coefficents and the fit covarience matrix
             #Coefficents are the slope and offset of the esat region
             esat_coeff, ccov = np.polyfit(esat_v, esat_i, 1, cov=True)
             esat_fit = esat_coeff[0]*voltage + esat_coeff[1]
             
             cerr = np.sqrt(np.diag(ccov))
             
             

             #Calculate indices based on the esat range
             exp_a = fixRange(voltage, exp_range[0])
             exp_b = fixRange(voltage, exp_range[1])
             
             #This indicates a bad curve that will not fit
             if exp_a > exp_b:
                 raise ValueError
             
             
             #Pull those parts of the curve
             exp_v = voltage[exp_a:exp_b]
             exp_i = current[exp_a:exp_b]
             
             
             #Fit the exponential region
             #three coefficients correspond to the amplitude, offset, and temperature
             popt, pcov = curve_fit(expFcn, exp_v, exp_i)
             A = popt[0]
             B = popt[1]
             kTe = popt[2]
             exp_fit = expFcn(voltage, A, B, kTe)
             
             
             perr = np.sqrt(np.diag(pcov))
             

             if plots:
                  fig_fit, ax_fit = plt.subplots(figsize = [8,4])
                  fig_fit.suptitle("Fitted Langmuir Curve")
                  ax_fit.set_xlabel("Sweep Voltage (V)")
                  ax_fit.set_ylabel("Electron Current (A)")
                  #Shade actual fitted regions
                  ax_fit.axvspan(voltage[exp_a], voltage[exp_b], alpha=0.15, color='red')
                  ax_fit.axvspan(voltage[esat_a], voltage[esat_b], alpha=0.25, color='green')
                  #Plot data, fix plot axes to match
                  ax_fit.plot(voltage, current, 'k', label='Data')
                  ax_fit.set_xlim(np.min(voltage), np.max(voltage))
                  ax_fit.set_ylim(np.min(current), np.max(current))
                  #Plot the fits
                  ax_fit.plot(voltage, esat_fit, 'g', linewidth=3, label='E. Sat. Fit')
                  ax_fit.plot(voltage, exp_fit, 'r', linewidth=3, label='Exp. Fit')
                  
                  ax_fit.legend()
                  
             #Calculate the plasma potential as the intersection of the exp and esat fits
             vpp_ind = np.argmin(np.abs(exp_fit - esat_fit))
             vpp = voltage[vpp_ind]
             #The actual error estimate for this quantity involes error propagation
             #through a difficult equation:
             #A e^{-x/kTe} + B - (C*x + D) = 0, solved for x
             
             
             #Calculate outputs
             esat = current[vpp_ind]
             
             vthe = 4.19e7*np.sqrt(kTe) #NRL formulay -> cm/s
             density = esat/(area*1.6e-19*vthe)
             
             #TODO This dummy error variable is a stand-in until someone 
             #can incorperate proper error propagation into this routine
             #There will be five error variables, corresponding to:
             #1,2 : slope and offset of esat
             #3,4: amplitude and offset of expoenntial region
             #5: temperature
             ebars = np.array([cerr[0], cerr[1], perr[0], perr[1], perr[2]])
             
             
             
             #Negative kTe indicates a bad fit
             if kTe < 0:
                 raise ValueError
             
             if verbose:
                  print("Vpp = " + str(vpp) + " V")
                  print("Te =  " + str(kTe) + " eV")
                  print("I_esat =  " + str(esat) + " A")
                  print("vthe =  " + str(vthe) + " cm/s")
                  print("n_e =  " + str(density) + " cm^-3")
                  
                  
         except(TypeError, RuntimeError, ValueError, np.linalg.LinAlgError):
             esat_fit = fail_val
             exp_fit = fail_val
             vpp = fail_val
             kTe = fail_val  
             esat = fail_val
             vthe = fail_val
             density=fail_val
             ebars = np.ones([5])*fail_val
             
           
     #Return
     if return_fits:
          return esat_fit, exp_fit, vpp, kTe, esat, vthe, density, ebars
     else:
          return vpp, kTe, esat, vthe, density, ebars
          
     



def vsweepLangmuirRawToFull(src, ndest, tdest, 
                  grid=True,
                  verbose=False, plots=False, debug=False,
                  grid_precision=0.1, strict_grid=False, strict_axes = False):
    """ Fits sweept Langmuir probe data and creates two full save files
    containing the calculated density and temperature

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing swept Langmuir probe data
            There should be two channels: the first being the Langmuir current
            and the second being the ramp voltage.
            
       ndest: hdfPath object
            Path string to location processed density data is written out
            
       tdest: hdfPath object
            Path string to location processed temperature data is written out
  
        grid: Boolean
            If grid is true, output will be written in cartesian grid array
            format, eg. [nti, nx, ny, nz, nreps, nchan]. Otherwise, output will
            be in [nshots, nti, nchan] format

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
        req_keys = ['area', 'resistor',
                    'gain', 'atten', 'ramp_gain', 'ramp_atten',
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z']

        if  'pos' in srcgrp:
            pos = srcgrp['pos'][:] #Read the entire array in
     
        else:
            #If no position information is given, a single explicit position
            #is required. 
            req_keys = req_keys + ['xpos', 'ypos', 'zpos']
            grid = False #Can't grid data if there's no pos array!

        #Process the required keys, throwing an error if any cannot be found
        csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        

        #Extract the shape of the source data
        nshots, nti, nchan = srcgrp['data'].shape
        
        #If requested by keyword, apply gridding
        if grid:
           shotlist, xaxis, yaxis, zaxis, nx, ny, nz, nreps, nshots = postools.grid(
                     pos, attrs, strict_axes=strict_axes, 
                     strict_grid=strict_grid, grid_precision=grid_precision, 
                     invert=True)

        if verbose:
            print("Opening destination HDF files")
        
        #Create the destination file directory if necessary
        #hdftools.requireDirs(ndest.file)
        #hdftools.requireDirs(tdest.file)

        #Open the destination file
        #This exists WITHIN the open statement for the source file, so the
        #source file is open at the same time.
        

        
       #Check if the output files exist already: if so, delete them
        if os.path.exists(ndest.file):
           os.remove(ndest.file)
        if os.path.exists(tdest.file):
           os.remove(tdest.file)
        

        with h5py.File(ndest.file, 'a') as ndf:
            with h5py.File(tdest.file, 'a') as tdf:
            
                #Throw an error if this group already exists
                if ndest.group != '/' and ndest.group in ndf.keys():
                     raise hdftools.hdfGroupExists(ndest)
                if tdest.group != '/' and tdest.group in tdf.keys():
                     raise hdftools.hdfGroupExists(tdest)
            
                ndestgrp = ndf.require_group(ndest.group)
                tdestgrp = tdf.require_group(tdest.group)
                
                
                grps = [ndestgrp, tdestgrp]
                
                for grp in grps:
                   hdftools.copyAttrs(srcgrp, grp)
                   #Throw an error if this dataset already exists
                   if 'data' in grp.keys():
                        raise hdftools.hdfDatasetExists(str(grp) + ' -> ' + "'data'")

                #Determine the time vector and nti
                #Assume first shot is representative of the vramp
                vramp = srcgrp['data'][0,:, 1] 
                time = srcgrp['time'][:]
                peaktimes, start, end = find_sweeps(time, vramp, plots=plots)
                nti = len(peaktimes)
                time = peaktimes
                    
                #Create the dataset 'data' appropriate to whether or not output
                #data will be gridded
                if verbose:
                    print("Creating 'data' group in destination file")
                if grid:
                    ndestgrp.require_dataset('data', (nti, nx, ny, nz), np.float32, chunks=True, compression='gzip')
                    ndestgrp.require_dataset('error', (nti, nx, ny, nz, 5), np.float32, chunks=True, compression='gzip')
                    
                    tdestgrp.require_dataset('data', (nti, nx, ny, nz), np.float32, chunks=True, compression='gzip')
                    tdestgrp.require_dataset('error', (nti, nx, ny, nz, 5), np.float32, chunks=True, compression='gzip')
                else:
                    ndestgrp.require_dataset('data', (nti,), np.float32, chunks=True, compression='gzip')
                    ndestgrp.require_dataset('error', (nti,5), np.float32, chunks=True, compression='gzip')
                    
                    tdestgrp.require_dataset('data', (nti,), np.float32, chunks=True, compression='gzip')
                    tdestgrp.require_dataset('error', (nti,5), np.float32, chunks=True, compression='gzip')


                
                probe_gain = float(attrs['gain'][0])
                probe_atten = float(attrs['atten'][0])
                ramp_gain = float(attrs['ramp_gain'][0])
                ramp_atten = float(attrs['ramp_atten'][0])
                
                resistor = float(attrs['resistor'][0]) #Ohms
                area = (attrs['area'][0]*u.Unit(attrs['area'][1])).to(u.cm ** 2).value
                
                probe_calib = np.power(10, probe_atten/20.0)/probe_gain
                ramp_calib = np.power(10, ramp_atten/20.0)/ramp_gain
                
                if grid:
                      #Initialize time-remaining printout
                      tr = util.timeRemaining(nx*ny*nz, reportevery=20)
                      for xi in range(nx):
                           for yi in range(ny):
                                for zi in range(nz):
                                    i = zi + yi*nz + xi*nz*ny
                                    
                                    if verbose:
                                         tr.updateTimeRemaining(i)
                                    
                                    s = shotlist[xi,yi,zi, :]
                                         
                                    current = srcgrp['data'][s,:, 0]*probe_calib/resistor
                                    vramp = srcgrp['data'][s,:, 1]*ramp_calib
                                    
                                    #Average over shots
                                    current = np.mean(current, axis=0)
                                    vramp = np.mean(vramp, axis=0)
                                    
                                    for ti in range(nti):
                                         a = int(start[ti])
                                         b = int(end[ti])

                                         vpp, kTe, esat, vthe, density, error = vsweep_fit(vramp[a:b], current[a:b], 
                                                      esat_range=None, exp_range=None, plots=False, area=area)
                                         
                                         ndestgrp['data'][ti, xi, yi, zi] = density
                                         ndestgrp['error'][ti, xi, yi, zi, :] = error
                                         
                                         tdestgrp['data'][ti, xi, yi, zi] = kTe
                                         tdestgrp['error'][ti, xi, yi, zi, :] = error
        
                else:  #Not gridded
                      current = srcgrp['data'][:,:, 0]
                      current = np.mean(current, axis=0)*probe_calib/resistor
                      vramp = srcgrp['data'][:,:, 1]
                      vramp = np.mean(vramp, axis=0)*ramp_calib
                      
                      for ti in range(nti):
                          a = start[ti]
                          b = end[ti]
                          vpp, kTe, esat, vthe, density, error = vsweep_fit(vramp[a:b], current[a:b], 
                                                      esat_range=None, exp_range=None, plots=False, area=area)
                          
                          ndestgrp['data'][ti] = density
                          ndestgrp['error'][ti, :] = error
                          
                          tdestgrp['data'][ti] = kTe
                          tdestgrp['error'][ti, :] = error
         
            
            
                for grp in grps:
                     #Write the axes as required by the format of the data written  
                     if grid:
                          grp.require_dataset('pos', (nshots, 3), np.float32, chunks=True)[:] = srcgrp['pos'][0:nshots]
                          for k in srcgrp['pos'].attrs.keys():
                              grp['pos'].attrs[k] = srcgrp['pos'].attrs[k]
          
                          dimlabels = ['time', 'xaxis', 'yaxis', 'zaxis']
                          grp.require_dataset('xaxis', (nx,), np.float32, chunks=True)[:] = xaxis
                          grp['xaxis'].attrs['unit'] = attrs['motion_unit'][0]
                          
                          grp.require_dataset('yaxis', (ny,), np.float32, chunks=True)[:] = yaxis
                          grp['yaxis'].attrs['unit'] = attrs['motion_unit'][0]
                          
                          grp.require_dataset('zaxis', (nz,), np.float32, chunks=True)[:] = zaxis
                          grp['zaxis'].attrs['unit'] = attrs['motion_unit'][0]
                     else:
                          dimlabels = ['time']

                     grp.require_dataset('time', (nti,), np.float32, chunks=True)
                     grp['time'][:] = time
                     grp['time'].attrs['unit'] = srcgrp['time'].attrs['unit']
                      
                     grp['data'].attrs['unit'] = 'G'
                                 
            
                ndestgrp['data'].attrs['unit'] = 'cm^{-3}'
                tdestgrp['data'].attrs['unit'] = 'eV'
                 
                ndestgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
                tdestgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]

            if verbose:
                print("End of Sweept Langmuir routine!")
                
            return True

          
     
     
def isatRawToFull(src, dest, 
                  ti=1.0, mu=4.0,
                  tdiode_hdf=None, grid=False,
                  offset_range=(0,100), offset_rel_t0 = (False, False), 
                  verbose=False, debug = False,
                  grid_precision=0.1, strict_grid=False, strict_axes = False):
    """ Integrates isat Langmuir probe data, calibrates output using information about the probe.

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing data
            
        dest: hdfPath object
            Path string to location processed data should be written out
            
        ti: Ion temperature (eV). Default assumption is 1 eV, which is typical
        of the LAPD LaB6 plasma. Scaling is as 1/sqrt(Ti).
        
        
        mu: Ion mass number (m_i/m_p = mu). Default is 4.0, for Helium.

        tdiode_hdf:  hdfPath object
            Path to a raw hdf5 file containing tdiode data. If no HDF file is
            provided, no timing correction will be applied.
            
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
        req_keys = ['area', 'atten', 'gain','resistor',
                    'dir', 'pol', 
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z',
                    'dt']
       


        if  'pos' in srcgrp:
            pos = srcgrp['pos'][:] #Read the entire array in
        else:
            #If no position information is given, a single explicit position
            #is required. 
            req_keys = req_keys + ['probe_xpos', 'probe_ypos', 'probe_zpos']
            grid = False #Can't grid data if there's no pos array!
            
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

        
            
        if verbose:
            print("Opening destination HDF file")
        
        #Create the destination file directory if necessary
        hdftools.requireDirs(dest.file)

        #Open the destination file
        #This exists WITHIN the open statement for the source file, so the
        #source file is open at the same time.
        
        #remove files if they already exist
        if os.path.exists(dest.file):
            os.remove(dest.file)
            
            
        with h5py.File(dest.file, 'a') as df:
            
            #Throw an error if this group already exists
            if dest.group is not '/' and dest.group in df.keys():
                raise hdftools.hdfGroupExists(dest)
            
            destgrp = df.require_group(dest.group)

            #Copy over attributes
            hdftools.copyAttrs(srcgrp, destgrp)



            #Load the time vector
            t = srcgrp['time']
            
            #If tdiode_hdf is set, load the pre-processed tdiode data
            if tdiode_hdf is not None:
                if verbose:
                    print("Loading tdiode array from file.")
                with h5py.File(tdiode_hdf.file, 'r') as sf:
                    grp = sf[tdiode_hdf.group]
                    t0indarr = grp['t0indarr'][:]
                    goodshots = grp['goodshots'][:]
                    tdiode_attrs = hdftools.readAttrs(grp)
                    
                #If tdiode was digitized with a different dt, this correction
                #will be necessary
                dt_ratio = float(attrs['dt'][0])/float(tdiode_attrs['dt'][0])
                t0indarr = (t0indarr/dt_ratio).astype(np.int32)
                
                #We will remove up to max_t0shift indices from each array such that
                #the t0 indices all line up.
                min_t0ind = np.min(t0indarr[goodshots])
                max_t0shift = np.max(t0indarr[goodshots]) - min_t0ind
                #Compute new nti
                nti = nti - max_t0shift 

                t = t[0:nti] - t[min_t0ind]




            #Throw an error if this dataset already exists
            if 'data' in destgrp.keys():
                    raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
                    
            #Create the dataset 'data' appropriate to whether or not output
            #data will be gridded
            if verbose:
                print("Creating 'data' group in destination file")
            if grid:
                destgrp.require_dataset('data', (nti, nx, ny, nz, nreps), np.float32, chunks=(np.min([nti, 20000]),1,1,1,1), compression='gzip')
            else:
                destgrp.require_dataset('data', (nshots, nti), np.float32, chunks=(1, np.min([nti, 20000])), compression='gzip')
            
            
    
          
            dt = ( attrs['dt'][0]*u.Unit(attrs['dt'][1])).to(u.s).value
            
            resistor = float(attrs['resistor'][0]) #Ohms
            area = (attrs['area'][0]*u.Unit(attrs['area'][1])).to(u.m ** 2).value
                
                

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

                else:
                    #By default, read in the entire dataset
                    ta = None
                    tb = None
          
                if debug:
                    print("Data range: [" + str(ta) + "," + str(tb) + "]")
   
                    
                
                #Read in the data from the source file
                voltage = srcgrp['data'][i,ta:tb, 0]
                
     
                #Calculate density
                #Equation is 2 from this paper: 10.1119/1.2772282
                #This is valid for the regime Te~Ti, which is approx true in
                #LAPD
                density = 1.6e9*np.sqrt(mu)*voltage/(resistor*area) #cm^-3
             

                if grid:
                    #Get location to write this datapoint from the shotgridind
                    xi = shotgridind[i, 0]
                    yi = shotgridind[i, 1]
                    zi = shotgridind[i, 2]
                    repi = shotgridind[i, 3]
                    #Write data
                    try:
                        destgrp['data'][:, xi, yi, zi, repi] = density
          
                    except ValueError as e:
                        print("ERROR!")
                        print(destgrp['data'].shape)
                        print(voltage.shape)
                        print([xi, yi, zi, repi])
                        raise(e)
                else:
                    #Write data
                    destgrp['data'][i,:] = density
                               

            if verbose:
                print("Writing axes to destination file")
            

            if grid:
                #Add the other axes and things we'd like in this file
                destgrp.require_dataset('pos', (nshots, 3), np.float32, chunks=True)[:] = srcgrp['pos'][0:nshots]
                for k in srcgrp['pos'].attrs.keys():
                    destgrp['pos'].attrs[k] = srcgrp['pos'].attrs[k]
                    
                dimlabels = ['time', 'xaxis', 'yaxis', 'zaxis', 'reps']
                
                destgrp.require_dataset('xaxis', (nx,), np.float32, chunks=True)[:] = xaxis
                destgrp['xaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
                destgrp.require_dataset('yaxis', (ny,), np.float32, chunks=True)[:] = yaxis
                destgrp['yaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
                destgrp.require_dataset('zaxis', (nz,), np.float32, chunks=True)[:] = zaxis
                destgrp['zaxis'].attrs['unit'] = attrs['motion_unit'][0]
                
                destgrp.require_dataset('reps', (nreps,), np.int32, chunks=True)[:] = np.arange(nreps)
                destgrp['reps'].attrs['unit'] = ''

            else:
                dimlabels = ['shots', 'time']
                
                destgrp.require_dataset('shots', (nshots,), np.int32, chunks=True)[:] = srcgrp['shots'][:]
                destgrp['shots'].attrs['unit'] = srcgrp['shots'].attrs['unit']
            
            
            destgrp.require_dataset('time', (nti,), np.float32, chunks=True)
            destgrp['time'][:] = t
            destgrp['time'].attrs['unit'] = srcgrp['time'].attrs['unit']


            destgrp['data'].attrs['unit'] = 'cm^{-3}'
                 
            destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            
            if verbose:
                print("End of isat Langmuir routine!")
                
            return True    

     
     
     
     
     
if __name__ == "__main__":
     
     """
     f=  hdftools.hdfPath(os.path.join("/Volumes", "PVH_DATA", "LAPD_Mar2018", "RAW","run104_JanusBaO_raw.hdf5"))
     ndest=  hdftools.hdfPath(os.path.join("/Volumes", "PVH_DATA", "LAPD_Mar2018", "FULL","run104_JanusBaO_density.hdf5"))
     tdest=  hdftools.hdfPath(os.path.join("/Volumes", "PVH_DATA", "LAPD_Mar2018", "FULL","run104_JanusBaO_temperature.hdf5"))
     
     #f=  hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "RAW","run104_JanusBaO_raw.hdf5"))
     #ndest=  hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL","run104_JanusBaO_density.hdf5"))
     #tdest=  hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL","run104_JanusBaO_temperature.hdf5"))
     resistor = 2.2
     
     vsweepLangmuirRawToFull(f, ndest, tdest, verbose=True, grid=True)
     
   
     
     """
     
     #f=  hdftools.hdfPath(os.path.join("/Volumes", "PVH_DATA", "LAPD_Mar2018", "RAW","run104_JanusLaB6_raw.hdf5"))
     f=  hdftools.hdfPath(os.path.join("F:",  "LAPD_Mar2018", "RAW","run104_JanusLaB6_raw.hdf5"))
     
     
    #isatRawToFull(f, ndest, verbose=True, grid=True)
 
     with h5py.File(f.file, 'a') as sf:
          resistor = 2.2
          shot = 400
          reps = 5
          time = sf['time'][:]*8
          current = sf['data'][shot:shot+reps,:,0]/resistor #2.2 Ohm resistor, convert to A
          voltage = sf['data'][shot:shot+reps,:,1]*100 #Attenuation was x100
          
          current = np.mean(current, axis=0)
          voltage = np.mean(voltage, axis=0)
     
     peaktime, a, b = choose_sweep(time, voltage, 0, plots=True)
     vpp, kTe, esat = vsweep_fit(voltage[a:b], current[a:b],
                plots=True, verbose=True)
     
     area = .1e-2 #cm^2
     vthe = 4.19e7*np.sqrt(kTe) #NRL formulay -> cm/s
     density = esat/(area*1.6e-19*vthe)
     
     print("Density = " + '{:.2e}'.format(density)+ " cm^-3")

    
          
     