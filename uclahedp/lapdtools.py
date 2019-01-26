#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
lapdtools.py: LAPD HDF5 analysis programs
 --> readRunPobe(run, probe, data_dir, dest, verbose=False)
     Reads LAPD HDF5 dataset into a new HDF5 file w/ metadata,
     calls lapdReadHDF with options chosen from metadata csvs.
        
--> lapdReadHDF(src=None, dest=None, channel_arr = None, 
         controls = None, verbose=False )
    Wrapper on bapsflib file.read_data which reads data from an LAPD HDF5 file
    and writes it out to an output HDF5 file.
"""

import numpy as np
from astropy import units as u
import time
import os

import hedpConstants as c
import csvtools as csvtools
import hdftools
import util
#bapsflib is available from pip. Run command 'pip bapsflib' from terminal to install it
from bapsflib import lapd

import h5py


def readRunProbe( run, probe, data_dir, dest, verbose=False):
    """ Retreives the appropriate metadata for a run and probe in a given data
    directory, then reads in the data using the bapsflib module and saves
    it in a new hdf5 file.
    
    Parameters
    ----------
        run: int
            Run number
        
        probe: str
            Probe name
            
        data_dir: str (path)
            Path to the data directory. Directory should have /HDF folder with 
            hdf files and /METADATA folder with CSV metadata files
    

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
    csv_dir = os.path.join(data_dir, c.metadata_dir)
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
    src =  os.path.join(data_dir, c.hdf_dir, attrs['datafile'][0] +  '.hdf5')
    
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
    
    #Check for the existance of a couple keys that define the existance of
    #an LAPD motion list. If they don't exist, the probe is assumed to be
    #stationary
    motion_keys = ['motion_controller', 'motion_receptacle']    
    if csvtools.missingKeys(attrs, motion_keys, fatal_error = False):
        print("Some motion keys not found: positon data will not be read out!")
        motion = False
    else:
        motion = True  
        
    #If motion=True, create the controls list for bapsflib
    # control = array of tuples of form (motion control, receptacle)
    # eg. controls = [('6K Compumotor', receptacle)]
    # note that 'receptacle' here is the receptacle NUMBER, 1 - indexed!)
    controls = None #Default value
    if motion:
        motion_controller = attrs['motion_controller'][0]
        motion_receptacle = attrs['motion_receptacle'][0]
        #If motion_controller isn't in this list, lapdReadHDF can't handle it
        if motion_controller in ['6K Compumotor', 'NI_XZ', 'NI_XYZ']:
            controls = [(motion_controller, motion_receptacle)]


    #Create the destination file and copy over attributes
    with h5py.File(dest.file, "a") as df:

        #Throw an error if this group already exists
        #This is to prevent memory allocation errors from overwriting
        if dest.group is not '/' and dest.group in df.keys():
            raise hdftools.hdfGroupExists(dest)
            
        #Create the group
        grp = df[dest.group]

        #Write the attrs dictioanry into attributes of the new data group
        hdftools.writeAttrs(attrs, grp)

    #Call lapdReadHDF to actually copy over the data
    lapdReadHDF(src=src, dest = dest, channel_arr = channel_arr, controls = controls, verbose=verbose)
    


def lapdReadHDF(src=None, dest=None, channel_arr = None, controls = None, verbose=False ):
    """ 
    Wrapper on the LAPD's bapsflib package that choses parametesr based on
    the metadata dictonary, manually handles several drives that aren't included
    in the bapsflib software, and writes out data in the standard HEDP
    HDF format
    
    Parameters
    ----------
        src: hdfPath object
            Path string to source hdf file (will be opend by bapsflib)
    

        dest: hdfPath object
            Path string to location data should be written out
            
        channel_arr: Array of tuples
            Array of channels formatted for input to the bapsflib package. 
            channel_arr = array of tuples of form:
                (digitizer, adc, board#, channel#)
            eg. channel_arr = ('SIS crate', 'SIS 3305', 2, 1)
            
            
        controls: Array of lists
            Array of motion controls formatted for input to the bapsflib 
            package. control = array of tuples of form:
            (motion control, receptacle)
            eg. controls = [('6K Compumotor', receptacle)]
            

        verbose: boolean
            Set this flag to true to enable print statements throughout the
            code, including a runtime-until-completion estimate during the
            data reading loop.

    Returns
    -------
       True, if execution is successful 
    """ 

    #This ugly initial code block extracts position information for drives
    #that are as-yet unsupported by bapsflib.
    #If the drive is unsupported, controls -> None so basflib doesn't look for it
    #bool(motion) keeps track of whether there WAS motion for later
    #Currently, this assumes only ONE XZ or XYZ drive is being used.
    #Not a problem for now, but just FYI
    if controls is not None:
        motion = True
        motion_attrs = {'unit':'cm'} #Units are assumed to be cm
        
        #If the controlling drive is the 6K Compumotor (XY drive)
        #This one is supported by bapsflib, so pos array will be
        #assembled later
        if controls[0][0] == '6K Compumotor':
            motion_format = 'fixed_rotation'
        
        #If the controlling drive is the NI_XZ drive...
        elif controls[0][0] == 'NI_XZ':
            controls = None
            motion_format = 'fixed_rotation' #Defines how to correct angles
            with h5py.File(src, "r") as sf:
                motion_group = sf['Raw data + config/NI_XZ']
                for item in motion_group.items():
                    #Find the group, which will be the motion group, and 
                    #copy over its attributes for future reference
                    if isinstance(item[1], h5py.Group):
                        config_name = str(item[0])
                        for k in motion_group.attrs.keys():
                            motion_attrs[k] = motion_group.attrs[k]
                        for k in motion_group[config_name].attrs.keys():
                            motion_attrs[k] = motion_group[config_name].attrs[k]
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
        elif controls[0][0] == 'NI_XYZ':
            motion_format = 'fixed_rotation'
            controls = None
            with h5py.File(src, "r") as sf:
                motion_group = sf['Raw data + config/NI_XYZ']
                for item in motion_group.items():
                    #Find the group, which will be the motion group, and 
                    #copy over its attributes for future reference
                    if isinstance(item[1], h5py.Group):
                        config_name = str(item[0])
                        for k in motion_group.attrs.keys():
                            motion_attrs[k] = motion_group.attrs[k]
                        for k in motion_group[config_name].attrs.keys():
                            motion_attrs[k] = motion_group[config_name].attrs[k]
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
   
    else:
        motion = False
        
        
    
    #Extract the number of shots from the HDF file
    #Currently seems to be no good way of doing this with bapsflib, at least
    #without loading the entire dataset into memory first... (Issue #24)
    with h5py.File(src, "r") as sf:
                run_seq = sf['Raw data + config/Data run sequence/Data run sequence']
                nshots = np.max(run_seq['Shot number'])
                
                
    #Determine the number of channels from the channel array
    nchan = len(channel_arr)

    #Load the source file as a bapsflib File object.
    sf = lapd.File(src, silent=True)

    #Read the first shot and channel to get the nti
    #We don't add controls here, because it causes an error in bapsflib 
    #(Issue #25)
    data = sf.read_data(channel_arr[0][2], channel_arr[0][3], 
                            digitizer =channel_arr[0][0], 
                            adc = channel_arr[0][1],
                            silent=True, shotnum=1)
   
    #Data['signal'] = (nshots, nti )
    nti = data['signal'].shape[1]

    #Read in digitizer parameters
    clock_rate = data.info['clock rate'].to(u.Hz)
    dt =  (  1.0 / clock_rate  ).to(u.s)
    
    
    #If controls is None, read the position information from the datafile
    if controls is not None:
        #Currently, bapsflib requires us to read in the entire dataset in order to
        #get the pos array. This isn't ideal, since it all goes into memory
        #We'll delete it asap and then chunk the rest...
        #Reading it within this if statement so we can avoid it unless
        #absolutely necessary
        data = sf.read_data(channel_arr[0][2], channel_arr[0][3], 
                            digitizer =channel_arr[0][0], 
                            adc = channel_arr[0][1],
                            add_controls=controls,
                            silent=True)
        pos =  data['xyz']
        del(data)
        
    



    #Open the destination file
    with h5py.File(dest.file, "a") as df:

        grp = df.require_group(dest.group)
        
        #Initialize the outpu data array
        if 'data' in grp.keys():
            raise hdftools.hdfDatasetExists(str(dest) + ' -> ' + "'data'")
        else:
            #Create the dataset + associated attributes
            grp.require_dataset("data", (nshots, nti, nchan), np.float32, 
                                chunks=(1, 20000, 1), compression='gzip')
            grp['data'].attrs['unit'] = 'V'
            grp.attrs['dt'] = [s.encode('utf-8') for s 
                     in [str(dt.value), str(dt.unit)] ]
            dimlabels = ['shots', 'time', 'chan']
            grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
            
            
        
        
        #Initialize some variables to use in time-remaining printout
        tstart  = time.time()
        tperstep = []
        nstepsreport = int(nshots*nchan/25.0)

        #Loop through the channels and shots, reading one-by-one into the
        #output dataset
        for chan in range(nchan):
            channel = channel_arr[chan]
            if verbose:
                print("Reading channel: " + str(chan+1) + '/' + str(nchan))
            

            for shot in range(nshots):
                #Update the time-per-shot time estimator
                #Over time this should give more accurate run-time
                #predictions
                nowtime = time.time()
                tperstep.append(nowtime-tstart)
                tstart = nowtime
                i = nshots*chan + shot
                if verbose and i % nstepsreport == 0 and i > 0:
                    tremain = np.mean(tperstep)*(nshots*nchan - i)
                    print(str(i) + '/' + str(nshots*nchan) + ' samples read, ' + 
                          util.timeFormat(tremain) + ' remaining' )
                
                #Read the data through bapsflib
                data = sf.read_data(channel[2], channel[3], digitizer =channel[0],
                                    adc = channel[1], 
                                    silent=True, shotnum=shot+1)
                
                grp['data'][shot,:,chan] = data['signal']
    
                

        #If motion=True, there should be a pos array, so write it to the file
        if motion:
            grp.require_dataset('pos', (nshots, 3), np.float32)[:] = pos
            del(pos)
            #Write any motion attributes to the pos array too
            for k in motion_attrs:
                grp['pos'].attrs[k] = motion_attrs[k]
                grp.attrs['motion_format'] = motion_format
            del(motion_attrs)
            
        
        #Create the axes
        grp.require_dataset('shots', (nshots,), np.float32, chunks=True )[:] = np.arange(nshots)
        grp['shots'].attrs['unit'] = ''
        
        t = np.arange(nti)*dt
        grp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t.value
        grp['time'].attrs['unit'] =  str(t.unit)
        
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = np.arange(nchan)
        grp['chan'].attrs['unit'] = ''
        
        
    #Clear the LAPD HDF file from memory
    del(sf, data, t)
   
    return True
    



    

    
    
#Simple test program
if __name__ == "__main__":

    data_dir = os.path.join("F:", "/LAPD_Mar2018/")

    dest = hdftools.hdfPath( r"F:/LAPD_Mar2018/RAW/test_save.hdf5")

    print('reading')
    util.mem()
    tstart = util.timeTest()
    x =  readRunProbe(102, 'tdiode', data_dir, dest, verbose=True)
    util.timeTest(t0=tstart)
    util.mem()
    print('done')
    
