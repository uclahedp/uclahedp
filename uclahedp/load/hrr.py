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
import os

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util

import h5py


def hrrToRaw( run, probe, hdf_dir, csv_dir, dest, verbose=False):
    """ Retreives the appropriate metadata for a run and probe in a given data
    directory, then reads in the data from the HRR hdf5 output file.
    
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
       dest (Filepath to destination file)
    """ 

    #Create a dictionary of attributes from the entire directory of CSV
    #files that applies to this probe and run
    attrs = csvtools.getAllAttrs(csv_dir, run, probe)
  
    #Check that some required keys are present, throw a fatal error if not
    req_keys = ['datafile']
    csvtools.missingKeys(attrs, req_keys, fatal_error=True)


    print(hdf_dir)
    print(attrs['datafile'])
    print(attrs['datafile'][0] +  '.hdf5')
    #TODO: Should this file take a data_dir and determine the filename
    #automatically, or should a source hdf file be given, leaving the program
    #that calls this one to determine the HDF file name?
    src =  os.path.join(hdf_dir,  attrs['datafile'][0] +  '.hdf5')
    
 
    #Create an array of channels
    #channel_arr = tuples of form (resource number, channel number)
    #Indexd from 1, to match load/LAPD.py
    channel_arr = []
    nchan = 1
    while True:
        digistr = 'resource' + str(int(nchan))
        chanstr = 'chan' + str(int(nchan))
        if chanstr in attrs.keys() and digistr in attrs.keys():
            #Check to make sure channel has actual non-nan values
            if not np.isnan(attrs[digistr][0]) and not np.isnan(attrs[chanstr][0]):
                #Append the channel to the list to be extracted
                channel_arr.append( (attrs[digistr][0], attrs[chanstr][0]) )
            nchan = nchan + 1
        else:
            break
    

    #Create a dictionary of position channels
    #channel_arr = tuples of form (resource number, channel number)
    ax = ['x', 'y', 'z']
    pos_chan = {}
    nchan = 1
    for i in range(3):
        digistr = ax[i] + 'pos_resource'
        chanstr = ax[i] + 'pos_chan'
        if chanstr in attrs.keys() and digistr in attrs.keys():
            #Check to make sure channel has actual non-nan values
            if not np.isnan(attrs[digistr][0]) and not np.isnan(attrs[chanstr][0]):
                #Append the channel to the list to be extracted
                pos_chan[ax[i]] = ( attrs[digistr][0], attrs[chanstr][0])
            else:
                pos_chan[ax[i]] = None
        else:
                pos_chan[ax[i]] = None
        
    
    print(pos_chan)


    #Determine the number of channels from the channel array
    nchan = len(channel_arr)
    
        
    #Read some variables from the src file
    with h5py.File(src, 'r')  as sf:
        
        digi_name = 'RESOURCE ' + str(channel_arr[0][0])
        print(digi_name)
        digigrp = sf[digi_name]
        

        resource_type = digigrp.attrs['RESOURCE TYPE'].decode('utf-8')
        
        attrs['RESOURCE ALIAS'] = (digigrp.attrs['RESOURCE ALIAS'].decode('utf-8'),'')
        attrs['RESOURCE DESCRIPTION'] = (digigrp.attrs['RESOURCE DESCRIPTION'].decode('utf-8'),'')
        attrs['RESOURCE ID'] = (digigrp.attrs['RESOURCE ID'],'')
        attrs['RESOURCE MODEL'] = (digigrp.attrs['RESOURCE MODEL'].decode('utf-8'),'')
        attrs['RESOURCE TYPE'] = (resource_type,'')
        resource_unit = digigrp['CHANNEL 0']['UNITS'][0].decode('utf-8')
        
        attrs['motion_unit'] = ('mm', '')
        

        if resource_type == 'SCOPE':
            dataname = 'TRACE'
            nshots = digigrp['CHANNEL 0'][dataname].shape[0]
            nti = digigrp['CHANNEL 0'][dataname].shape[1]
            dt = digigrp['CHANNEL 0'][dataname].attrs['WAVEFORM DT'] * u.s
            
            attrs['dt'] = [str(dt.value), str(dt.unit) ]
            
            #attrs['dt'] = [s.encode('utf-8') for s 
            #     in [str(dt.value), str(dt.unit) ] ]
            
        elif resource_type == 'MOTOR BOARD':
            dataname = 'POSITION'
            nshots = digigrp['CHANNEL 0'][dataname].shape[0]
            nti = 1
        
        


    #Create the destination file
    with h5py.File(dest.file, "a") as df:

        #Create the dest group, throw error if it exists
        if dest.group is not '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
        grp = df[dest.group]
        
        #Initialize the output data array
        if 'data' in grp.keys():
            raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
            
        #Create the dataset + associated attributes
        grp.require_dataset("data", (nshots, nti, nchan), np.float32, 
                            chunks=(1, np.min([nti, 20000]), 1), 
                            compression='gzip')
        grp['data'].attrs['unit'] = resource_unit
        

        
        dimlabels = ['shots', 'time', 'chan']
        
        grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
        
    

        #Open the hdf5 file and copy the data over
        with h5py.File(src) as sf:
            
            
            #Initialize time-remaining printout
            tr = util.timeRemaining(nchan*nshots)
            
            #Loop through the channels and shots, reading one-by-one into the
            #output dataset
            for chan in range(nchan):
                digi_name = 'RESOURCE ' + str(channel_arr[chan][0])
                chan_name = 'CHANNEL ' + str(channel_arr[chan][1])
                if verbose:
                    print("Reading channel: " + str(chan+1) + '/' + str(nchan))
                for shot in range(nshots):
                    
                    if verbose:
                        tr.updateTimeRemaining(nshots*chan + shot)
                    #Read the data from the hdf5 file
                    grp['data'][shot,:,chan] = sf[digi_name][chan_name][dataname][shot, ...]
                    
                    
            
            if pos_chan['x'] is not None or pos_chan['y'] is not None or pos_chan['z'] is not None:
                
                grp.require_dataset('pos', (nshots, 3), np.float32)
                ax = ['x', 'y', 'z']
                
                unit_factor = (1.0*u.Unit(attrs['motion_unit'][0])).to(u.cm).value
                attrs['motion_unit'] = ('cm', '')

                for i, a in enumerate(ax):
                    if pos_chan[a] is not None:
                        resname = 'RESOURCE ' + str(pos_chan[a][0])
                        channame = 'CHANNEL ' + str(int(pos_chan[a][1]))
                        
                        
                        posdata = sf[resname][channame]['POSITION'][:]*unit_factor
                        
                        #Handle the case where the multiple data points were
                        #taken at a position so npos!=nshots
                        npos = posdata.size
                        if npos != nshots:
                            posdata = np.repeat(posdata, int(nshots/npos))

                        grp['pos'][:, i] = posdata
                    
                    else:
                        grp['pos'][:, i] = np.zeros(nshots)

        
        #Create the axes
        grp.require_dataset('shots', (nshots,), np.float32, chunks=True )[:] = np.arange(nshots)
        grp['shots'].attrs['unit'] = ''
        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ''
        
        if resource_type == 'SCOPE':
            t = np.arange(nti)*dt
            grp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t.value
            grp['time'].attrs['unit'] =  str(t.unit)
            
            
        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)
        
    return dest
    

#Simple test program
if __name__ == "__main__":

    #hdf_dir = os.path.join("F:", "LAPD_Mar2018", "HDF")
    #csv_dir = os.path.join("F:", "LAPD_Mar2018", "METADATA")
    #dest = hdftools.hdfPath( r"F:/LAPD_Mar2018/RAW/run102_PL11B_raw.hdf5")
    
    
    probe = 'tdiode'
    #probe = 'PLL_B1'
    run = 6
    
    
    hdf_dir = '/Volumes/PVH_DATA/2019BIERMANN/HDF/'
    csv_dir = '/Volumes/PVH_DATA/2019BIERMANN/METADATA/'
    dest = hdftools.hdfPath( '/Volumes/PVH_DATA/2019BIERMANN/RAW/' + 'run' + str(run) + '_' + probe + '_raw.hdf5')

    #Delete the output file if it already exists
    try:
        os.remove(dest.file)
    except FileNotFoundError:
        pass
    
    print('reading')
    util.mem()
    tstart = util.timeTest()
    x =  hrrToRaw(run, probe, hdf_dir, csv_dir, dest, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')

