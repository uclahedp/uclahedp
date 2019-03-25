#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
hrrtools.py: HEDP group High Rep Rate Program tools
 --> readRunPobe(run, probe, data_dir, dest, verbose=False)
     Reads an HDF5 file from the labview HRR program 

"""

import numpy as np
from astropy import units as u
import time
import os

from uclahedp import csvtools, hdftools, util

import h5py


def hrrToRaw( run, probe, hdf_dir, csv_dir, dest, verbose=False):
    """ Retreives the appropriate metadata for a run and probe in a given data
    directory, then reads in the data using the bapsflib module and saves
    it in a new hdf5 file.
    
    Parameters
    ----------
        run: int
            Run number
        
        probe: str
            Probe name
            
        hdf_dir: str (path)
            Path to the directory where HDF files are stored
            
        csv_dir: str(path)
            Path to the directory where metadata CSV's are stored
    

        dest: hdfPath object
            Path string to location data should be written out

        verbose: boolean
            Set this flag to true to enable print statements throughout the
            code, including a runtime-until-completion estimate during the
            data reading loop.

    Returns
    -------
       True, if execution is successful 
    """ 

    #Create a dictionary of attributes from the entire directory of CSV
    #files that applies to this probe and run
    attrs = csvtools.getAllAttrs(csv_dir, run, probe)
  
    #Check that some required keys are present, throw a fatal error if not
    req_keys = ['datafile', 'digitizer']
    csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        
    #Load some digitizer parameters we now know exist
    digitizer = attrs['digitizer'][0]
    digitizer = 'RESOURCE ' + str(digitizer)

    #TODO: Should this file take a data_dir and determine the filename
    #automatically, or should a source hdf file be given, leaving the program
    #that calls this one to determine the HDF file name?
    src =  os.path.join(hdf_dir,  attrs['datafile'][0] +  '.hdf5')
    
 
    #Create an array of channels (required input for bapsflib read_data)
    # channel_arr = array of tuples of form: (digitizer, adc, board#, channel#)
    # eg. channel_arr = ('SIS crate', 'SIS 3305', 2, 1)
    #Do this in a loop, so the number of channels is flexible
    #However, the number of 'brd' and 'chan' fields MUST match
    #AND, the keys must be of the format 'brd1', 'chan1', etc.
    channel_arr = []
    nchan = 1
    while True:
        chanstr = 'chan' + str(int(nchan))
        if chanstr in attrs.keys():
            #Check to make sure channel has actual non-nan values
            if not np.isnan(attrs[chanstr][0]):
                #Append the channel to the list to be extracted
                channel_arr.append( (digitizer, adc, attrs[brdstr][0], attrs[chanstr][0]) )
            nchan = nchan + 1
        else:
            break
    #Determine the number of channels from the channel array
    nchan = len(channel_arr)
        
    #Read some variables from the src file
    with h5py.File(src, silent=True)  as sf:
        
        digigrp = sf[digitizer]
        
        attrs['RESOURCE_ALIAS'] = digigrp.attrs['RESOURCE_ALIAS']
        
        #TODO: Read the rest of these attributes here
        dt = digigrp['CHANNEL 0']['TRACE']['WAVEFORM DT']
        nshots = digigrp['CHANNEL 0']['TRACE'].shape[0]
        nti = digigrp['CHANNEL 0']['TRACE'].shape[1]
        digi_unit = digigrp['CHANNEL 0']['UNITS'][0]



    #Create the destination file
    with h5py.File(dest.file, "a") as df:

        #Create the dest group, throw error if it exists
        if dest.group is not '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
        grp = df[dest.group]
        
        #Initialize the output data array
        if 'data' in grp.keys():
            raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
        
        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)

        #Open the hdf5 file and copy the data over
        with h5py.File(src, silent=True) as sf:
            
            digigrp = sf[digitizer]
            
            #Create the dataset + associated attributes
            grp.require_dataset("data", (nshots, nti, nchan), np.float32, 
                                chunks=(1, np.min([nti, 20000]), 1), 
                                compression='gzip')
            grp['data'].attrs['unit'] = digi_unit
            
            grp.attrs['dt'] = [s.encode('utf-8') for s 
                     in [str(dt.value), str(dt.unit)] ]
            
            dimlabels = ['shots', 'time', 'chan']
            
            grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
                

            #Initialize time-remaining printout
            tr = util.timeRemaining(nchan*nshots)
            
            #Loop through the channels and shots, reading one-by-one into the
            #output dataset
            for chan in range(nchan):
                channel = channel_arr[chan]
                if verbose:
                    print("Reading channel: " + str(chan+1) + '/' + str(nchan))

                for shot in range(nshots):
                    
                    if verbose:
                        tr.updateTimeRemaining(nshots*chan + shot)
                    
                    #Read the data from the hdf5 file
                    grp['data'][shot,:,chan] = digigrp['CHANNEL' + str(chan)]['TRACE'][shot, :]

        #Create the axes
        grp.require_dataset('shots', (nshots,), np.float32, chunks=True )[:] = np.arange(nshots)
        grp['shots'].attrs['unit'] = ''
        
        t = np.arange(nti)*dt
        grp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t.value
        grp['time'].attrs['unit'] =  str(t.unit)
        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ''

    del(sf, data, t)
   
    return dest
    

#Simple test program
if __name__ == "__main__":

    #hdf_dir = os.path.join("F:", "LAPD_Mar2018", "HDF")
    #csv_dir = os.path.join("F:", "LAPD_Mar2018", "METADATA")
    #dest = hdftools.hdfPath( r"F:/LAPD_Mar2018/RAW/run102_PL11B_raw.hdf5")
    
    hdf_dir = '/Volumes/PVH_DATA/2019BIERMANN/HDF/'
    csv_dir = '/Volumes/PVH_DATA/2019BIERMANN/METADATA/'
    dest = hdftools.hdfPath( '/Volumes/PVH_DATA/2019BIERMANN/RAW/')

    print('reading')
    util.mem()
    tstart = util.timeTest()
    x =  lapdToRaw(102, 'PLL_B1', hdf_dir, csv_dir, dest, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')

