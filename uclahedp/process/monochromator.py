# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 15:14:46 2019

@author: jessi

--> monochromatorRawToFull(src, dest, tdiode_hdf=None, grid=False, verbose=False)
    Takes in a source HDF5 file
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



def monochromatorRawToFull(src, dest, port=14, tdiode_hdf=None,
                  verbose=False, debug = False, vdist=False):
    """ 

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing bdot data
            
        dest: hdfPath object
            Path string to location processed bdot data should be written out

        tdiode_hdf:  hdfPath object
            Path to a raw hdf5 file containing tdiode data. If no HDF file is
            provided, no timing correction will be applied. 
            
        port: float
            port at which the probe is located
            
             
            
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
        req_keys = []
       
        
        #Process the required keys, throwing an error if any cannot be found
        csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        

        #Extract the shape of the source data
        nshots, nti, nchan = srcgrp['data'].shape

            
            
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
        
            destgrp.require_dataset('data', (nshots, nti), np.float32, chunks=(1, np.min([nti, 20000])), compression='gzip')
            
            #Load the time vector
            t = srcgrp['time']
            #If a timing diode is being applied, correct the time vector here.
            if tdiode_hdf is not None:
                t = t[0:nti] - t[min_t0ind]
    
         

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
                signal = - np.squeeze(srcgrp['data'][i,ta:tb])
                
                if vdist:
                    destgrp['data'][i,:] = np.flip(signal)
                else:
                    destgrp['data'][i,:] = signal
                    
            destgrp['data'].attrs['unit'] = ''
                       
            destgrp.require_dataset('shots', (nshots,), np.int32, chunks=True)[:] = srcgrp['shots'][:]
            destgrp['shots'].attrs['unit'] = srcgrp['shots'].attrs['unit']
            

            if vdist:
                dist = (port-13)*.325    
                v = np.where(t!=0, dist/t,0)
                dimlabels = ['shots', 'velocity']
                destgrp.require_dataset('velocity', (nti,), np.float32, chunks=True)[:] = np.flip(v)
                destgrp['velocity'].attrs['unit'] = 'm/s'
            else:
                dimlabels = ['shots', 'time']
                destgrp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t
                destgrp['time'].attrs['unit'] = 's'
            
            destgrp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            

            if verbose:
                print("End of Monochromator routine!")
                
            return True
    

if __name__ == "__main__":
     
     csvfile = os.path.join("F:","LAPD_Mar2018","Bdot Calibration Data", "LAPD7.csv")
     #csvfile = os.path.join("/Volumes", "PVH_DATA","LAPD_Jan2019","bdot_calibrations", "LAPD_C6.csv")
     #csvfile = os.path.join("/Volumes", "PVH_DATA","LAPD_Mar2018","Bdot Calibration Data", "LAPD7.csv")
     
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
