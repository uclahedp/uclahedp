#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
lapd.py: LAPD HDF5 analysis programs
 --> readRunPobe(run, probe, data_dir, dest, verbose=False)
     Reads LAPD HDF5 dataset into a new HDF5 file w/ metadata,
     calls lapdReadHDF with options chosen from metadata csvs.

"""

import numpy as np
from astropy import units as u
import time
import os

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util

# bapsflib is available from pip. Run command 'pip bapsflib' from terminal to
# install it
from bapsflib import lapd as bapsf_lapd

import h5py


def lapdToRaw( run, probe, hdf_dir, csv_dir, dest, verbose=False):
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
    req_keys = ['datafile', 'digitizer', 'adc']
    csvtools.missingKeys(attrs, req_keys, fatal_error=True)
        
    #Load some digitizer parameters we now know exist
    digitizer = attrs['digitizer'][0]
    adc = attrs['adc'][0]
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
        brdstr = 'brd' + str(int(nchan))
        chanstr = 'chan' + str(int(nchan))
        if brdstr in attrs.keys() and chanstr in attrs.keys():
            #Check to make sure channel has actual non-nan values
            if not np.isnan(attrs[brdstr][0])  and not np.isnan(attrs[chanstr][0]):
                #Append the channel to the list to be extracted
                channel_arr.append( (digitizer, adc, attrs[brdstr][0], attrs[chanstr][0]) )
            nchan = nchan + 1
        else:
            break
    #Determine the number of channels from the channel array
    nchan = len(channel_arr)
        
    #Read some variables from the src file
    with bapsf_lapd.File(src, silent=True)  as sf:
        src_digitizers = sf.digitizers

        digi = src_digitizers['SIS crate'] #Assume this id the digitizer: it is the only one
        #Assume the adc, nti, etc. are all the same on all the channels.
        
        #This line assumes that only one configuration is being used
        #This is usually the case: if it is not, changes need to be made
        daq_config = digi.active_configs[0]
        
        name, info = digi.construct_dataset_name(channel_arr[0][2], 
                                                 channel_arr[0][3], 
                                                 adc=channel_arr[0][1],
                                                 config_name = daq_config,
                                                 return_info=True)
        #Read out some digitizer parameters
        nshots = info['nshotnum']
        nti = info['nt']
        clock_rate = info['clock rate'].to(u.Hz)
        dt =  (  1.0 / clock_rate  ).to(u.s)

    
    #Check if keys are provided to specify a motion list
    # control = array of tuples of form (motion control, receptacle)
    # eg. controls = [('6K Compumotor', receptacle)]
    # note that 'receptacle' here is the receptacle NUMBER, 1 - indexed!)
    req_keys = ['motion_controller', 'motion_receptacle']    
    if csvtools.missingKeys(attrs, req_keys, fatal_error = False):
        print("Some motion keys not found: positon data will not be read out!")
        controls, pos = None, None
    else:
        motion_controller = attrs['motion_controller'][0]
        motion_receptacle = attrs['motion_receptacle'][0]
        controls = [(motion_controller, motion_receptacle)]
        
        #Check to see if the motion controller reported actually exists in the
        #hdf file. If not, assume the probe was stationary (motion=None)
        #If motion_controller isn't in this list, lapdReadHDF can't handle it
        #Check if the motion controller provided is supported by the code and
        if motion_controller in ['6K Compumotor', 'NI_XZ', 'NI_XYZ']:
            pos, attrs = readPosArray(src, controls, attrs)
        else:
            controls, pos = None, None
    
    #Create the destination file directory if necessary
    hdftools.requireDirs(dest.file)
    #Create the destination file
    with h5py.File(dest.file, "a") as df:

        #Create the dest group, throw error if it exists
        if dest.group is not '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
        grp = df[dest.group]
        
        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)

        #Open the LAPD file and copy the data over
        with bapsf_lapd.File(src, silent=True) as sf:
            
            #Initialize the output data array
            if 'data' in grp.keys():
                raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")

            #Create the dataset + associated attributes
            grp.require_dataset("data", (nshots, nti, nchan), np.float32, 
                                chunks=(1, np.min([nti, 20000]), 1), 
                                compression='gzip')
            grp['data'].attrs['unit'] = ('V', '')
            
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
                    
                    #Read the data through bapsflib
                    data = sf.read_data(channel[2], channel[3], digitizer =channel[0],
                                        adc = channel[1], config_name = daq_config, 
                                        silent=True, shotnum=shot+1)
                    
                    grp['data'][shot,:,chan] = data['signal']
    
                

        #If applicable, write the pos array to file
        if pos is not None:
            grp.require_dataset('pos', (nshots, 3), np.float32)[:] = pos
            del pos
            
        
        #Create the axes
        grp.require_dataset('shots', (nshots,), np.float32, chunks=True )[:] = np.arange(nshots)
        grp['shots'].attrs['unit'] = ('','')
        
        t = np.arange(nti)*dt
        grp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t.value
        grp['time'].attrs['unit'] =  (str(t.unit), '')
        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ('', '')
        
        
    #Clear the LAPD HDF file from memory
    del(sf, data, t)
   
    return dest
    

def readPosArray(src, controls, motion_attrs):
     
    #This ugly initial code block extracts position information for drives
    #that are as-yet unsupported by bapsflib.
    #bool(motion) keeps track of whether there WAS motion for later
    #Currently, this assumes only ONE XZ or XYZ drive is being used.
    #Not a problem for now, but just FYI
    motion_attrs['motion_unit'] = (str("cm"), '')
    
    print("Reading position array: " + str(controls[0][0]))
    
    control = controls[0][0]
    receptacle = controls[0][1]

    #For the 6K, we can check to see if there is actually one present
    if control == '6K Compumotor' :
        motion_attrs['motion_format'] = ('fixed_pivot', '')  
        with bapsf_lapd.File(src, silent=True)  as sf: 
            
            if control in sf.controls.keys():
                src_controls = sf.controls
                pos = sf.read_controls(controls)['xyz']
                
                #Extract some other parameters
                #This is currently ALL the motion lists in the file
                temp_list = sf.controls['6K Compumotor'].configs[receptacle]['motion lists']
                dset= sf[sf.controls['6K Compumotor'].configs[receptacle]['dset paths'][0]]
                motion_list = ( dset['Motion list'][0] ).decode("utf-8")
                
                nx = temp_list[motion_list]['npoints'][0]
                ny = temp_list[motion_list]['npoints'][1]
                nz = temp_list[motion_list]['npoints'][2]
                
                dx = temp_list[motion_list]['delta'][0]
                dy = temp_list[motion_list]['delta'][1]
                dz = temp_list[motion_list]['delta'][2]
                
                x0 = temp_list[motion_list]['center'][0]
                y0 = temp_list[motion_list]['center'][1]
                z0 = temp_list[motion_list]['center'][2]
                

    #If the controlling drive is the NI_XZ drive...
    elif control == 'NI_XZ':
        motion_attrs['motion_format'] = ('fixed_pivot', '')   #Defines how to correct angles
        with h5py.File(src, "r") as sf:
            motion_group = sf['Raw data + config/NI_XZ']
            for item in motion_group.items():
                if isinstance(item[1], h5py.Group):
                    motion_list = str(item[0])
                    mlpath = 'Raw data + config/NI_XZ/' + motion_list

                    nx, ny, nz = sf[mlpath].attrs['Nx'], 1.0, sf[mlpath].attrs['Nz']
                    dx, dy, dz = sf[mlpath].attrs['dx'], 0.0, sf[mlpath].attrs['dz']
                    x0, y0, z0 = sf[mlpath].attrs['x0'], 0.0, sf[mlpath].attrs['z0']
                     
    
            runtimelist = sf['Raw data + config/NI_XZ/Run time list']
            xpos =   runtimelist['x']
            zpos =   runtimelist['z']
            nshots = len(xpos)
            ypos = np.zeros([nshots])
            pos = np.zeros([nshots, 3])
            pos[:,0] = xpos
            pos[:,1] = ypos
            pos[:,2] = zpos
            del(xpos, ypos, zpos)

            
    #If the controlling drive is the NI_XYZ drive (Krishna's)...
    elif control == 'NI_XYZ':
        motion_attrs['motion_format'] = ('fixed_pivot', '')  
        with h5py.File(src, "r") as sf:
            try:
                motion_group = sf['Raw data + config/NI_XYZ']
            except KeyError:
                return None
            
            for item in motion_group.items():
                if isinstance(item[1], h5py.Group):
                    motion_list = str(item[0])
                    mlpath = 'Raw data + config/NI_XYZ/' + motion_list

                    nx, ny, nz = sf[mlpath].attrs['Nx'], 1.0, sf[mlpath].attrs['Nz']
                    dx, dy, dz = sf[mlpath].attrs['dx'], 0.0, sf[mlpath].attrs['dz']
                    x0, y0, z0 = sf[mlpath].attrs['x0'], 0.0, sf[mlpath].attrs['z0']


            runtimelist = sf['Raw data + config/NI_XYZ/Run time list']
            xpos =   runtimelist['x']
            ypos = runtimelist['y']
            zpos =   runtimelist['z']
            nshots = len(xpos)
            pos = np.zeros([nshots, 3])
            pos[:,0] = xpos
            pos[:,1] = ypos
            pos[:,2] = zpos
            del(xpos, ypos, zpos)
            
    motion_attrs['motion_list'] = (motion_list, '')
    motion_attrs['nx'] = (nx, "")
    motion_attrs['ny'] = (ny, "")
    motion_attrs['nz'] = (nz, "")
    motion_attrs['dx'] = (dx, "")
    motion_attrs['dy'] = (dy, "")
    motion_attrs['dz'] = (dz, "")
    motion_attrs['x0'] = (x0, "")
    motion_attrs['y0'] = (y0, "")
    motion_attrs['z0'] = (z0, "")

    return pos, motion_attrs

    

    
    
#Simple test program
if __name__ == "__main__":
    
    exp = 'LAPD_Jan2019'
    probe = 'LAPD_C6'
    run = 30

    
    hdf_dir =  os.path.join("F:", exp, "HDF")
    csv_dir =  os.path.join("F:", exp, "METADATA")
    dest = hdftools.hdfPath( os.path.join("F:", exp, "RAW", 'run' + str(run) + '_' + probe + '_raw.hdf5'))

    hdf_dir = '/Volumes/PVH_DATA/' + exp + '/HDF/'
    csv_dir = '/Volumes/PVH_DATA/' + exp + '/METADATA/'
    dest = hdftools.hdfPath( '/Volumes/PVH_DATA/' +  exp +  '/RAW/run' + str(run) + '_' + probe + '_raw.hdf5')
    
    
    #Delete the output file if it already exists
    try:
        os.remove(dest.file)
    except FileNotFoundError:
        pass

    print('reading')
    util.mem()
    tstart = util.timeTest()
    x =  lapdToRaw(run, probe, hdf_dir, csv_dir, dest, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
    
