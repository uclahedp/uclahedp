#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program opens an LAPD HDF file (as well as metadata csvs)
and creates an ndfFile object containing the data. 

@author: peter
"""

#TODO currently 'NI_XZ' isn't valid drive

import numpy as np
from astropy import units as u

import os

import hedpConstants as c
import csvtools as csvtools
#bapsflib is available from pip. Run command 'pip papsflib' from terminal to install it
from bapsflib import lapd
from bapsflib._hdf import HDFMap



# channel_arr = array of tuples of form: (digitizer, adc, board#, channel#)
# eg. channel_arr = ('SIS crate', 'SIS 3305', 2, 1)
#
# control = array of tuples of form (motion control, receptacle)
# eg. controls = [('6K Compumotor', receptacle)]
# note that 'receptacle' here is the receptacle NUMBER, 1 - indexed!)
# if no control is specified, the measurement is assumed to be stationary
def bapsfReadHDF(src=None, dest=None, channel_arr = None, controls = None ):
    
    motion = not controls == None
    nchan = len(channel_arr)
    print('motion = ' + str(motion))

    f = lapd.File(src, silent=True)
    
    arr = [] # Array for temporarily storing the read data
    for i in range(nchan):
        channel = channel_arr[i]
        data = f.read_data(channel[2], channel[3],
                           digitizer =channel[0], adc = channel[1],
                           add_controls=controls, 
                           silent=True, )
        
        
        signal = data['signal'].T
        
        # Only need to do this for one channel, since we assume
        # that all channels coorespond to a single physical probe
        if i == 0:
            shotnum = data['shotnum']
            nshots = len(shotnum)
            nti = signal.shape[0]
            clock_rate = data.info['clock rate'].to(u.Hz)
            dt =  (  1.0 / clock_rate  ).to(u.s)
            t = np.arange(nti)*dt
            chan_axis = np.arange(nchan)*u.Unit('')
            shots_axis = np.arange(nshots)*u.Unit('')
            

         # Reshape the array according to the axes read in
        if motion and gridded:
            signal = np.reshape(signal, (nti, nx, ny, nz, nreps))
        else:
            signal = np.reshape(signal, (nti, nshots))
                
       
        # Temporarily store the array
        arr.append(signal)
        
    
    axes = []
    
    # Always include a time axis in the beginning
    if nti > 1:
        axes.append( {'label':'t', 'axis':t} )
    
    
    #Deal with motional gridded datasets
    if motion and gridded:
        output = np.zeros([nti, nx, ny, nz, nreps, nchan])

        if nx > 1:
            axes.append( {'label':'x', 'axis':xaxis} )
        if ny > 1:
            axes.append( {'label':'y', 'axis':yaxis} )
        if nz > 1:
            axes.append( {'label':'z', 'axis':zaxis} )
        if nreps > 1:
            axes.append( {'label':'reps', 'axis':rep_axis} )

    #Deal with motional non-gridded datasets
    #Deal with non-motional datasets (those with no motion device set  
    #(here both of these options have the same axes and shape requirements)
    else: 
        output = np.zeros([nti, nshots, nchan])
        if nshots > 1:
            axes.append( {'label':'shots', 'axis':shots_axis} )
      
        
        
        
    #Always include a channel axis at the end if there is more than one   
    if nchan > 1:
            axes.append( {'label':'channels', 'axis':chan_axis} )
            
            
    # Create an output data array with the proper form
    
    for i in range(nchan):
        output[..., i] = arr[i]
        
    # Eliminate trivial dimensions 
    output = np.squeeze(output)
    # Set units of data to volts
    output = output*u.V
    # Set datalabel
    data_label = 'Raw LAPD data'
    
   
    
    #TODO add the pos array to the ndf object as an option and save as that type here
    obj = sdfarr(data=output, axes=axes, data_label = data_label)
    obj.appendLog('Created by bapsftools.bapsfReadHDF from HDF file: ' + filepath)
    
    return obj
  
    


def readRunProbe( run, probe, exp_dir, dest):

    csv_dir = os.path.join(exp_dir, c.metadata_dir)

    run_level_attrs = getRunLevelAttrs(csv_dir, run)
    probe_level_attrs =  getProbeLevelAttrs(csv_dir, run, probe)
    
    attrs = {**run_level_attrs,  **probe_level_attrs}

    
    req_keys = ['datafile', 'digitizer', 'adc']
    motion_keys = ['motion_controller', 'motion_receptacle']
    
    missing_keys = []
    for k in req_keys:
        if not k in attrs.keys():
            missing_keys.append(k)
    if len(missing_keys) > 0:
        raise ValueError("Missing columns in csv files! The following keys were not found, and are required: " + str(missing_keys))
        
        
    missing_keys = []
    for k in motion_keys:
        if not k in attrs.keys():
            missing_keys.append(k)
    if len(missing_keys) > 0:
        print("Some motion keys not found: positon data will not be read out!")
        motion = False
    else:
        motion = True


    digitizer = attrs['digitizer'][0]
    adc = attrs['adc'][0]
    channel_arr = []
    nchan = 1
    while True:
        brdstr = 'brd' + str(int(nchan))
        chanstr = 'chan' + str(int(nchan))
        if brdstr in attrs.keys() and chanstr in attrs.keys():
            channel_arr.append( (digitizer, adc, attrs[brdstr][0], attrs[chanstr][0]) )
            nchan = nchan + 1
        else:
            break

        
    if motion:
        motion_controller = attrs['motion_controller'][0]
        motion_receptacle = attrs['motion_receptacle'][0]
        controls = [(motion_controller, motion_receptacle)]
    
    print(channel_arr)
    print(controls)
    
    datafile = attrs['datafile'][0]

    obj = bapsfReadHDF(src=datafile, dest = dest, channel_arr = channel_arr, controls = controls)

    return obj

    
    
    
if __name__ == "__main__":
    
    exp_dir = os.path.join("F:", "/LAPD_Mar2018/")


    src  = r"/F:/LAPD_Mar2018/HDF/peening074_apr09.hdf5" #Small file
    #sfile = r"/F:/LAPD_Mar2018/HDF/peening065_apr06.hdf5" #big file
    
    dest = r"/F:/LAPD_Mar2018/RAW/test_save.hdf5"

    print('reading')
    x =  readRunProbe(102, 'PL11B', exp_dir, dest)

    print('done')
