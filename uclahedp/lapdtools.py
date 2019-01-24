#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program opens an LAPD HDF file (as well as metadata csvs)

@author: peter
"""


import numpy as np
from astropy import units as u

import os

import hedpConstants as c
import csvtools as csvtools
import hdftools
#bapsflib is available from pip. Run command 'pip papsflib' from terminal to install it
from bapsflib import lapd

import h5py


# channel_arr = array of tuples of form: (digitizer, adc, board#, channel#)
# eg. channel_arr = ('SIS crate', 'SIS 3305', 2, 1)
#
# control = array of tuples of form (motion control, receptacle)
# eg. controls = [('6K Compumotor', receptacle)]
# note that 'receptacle' here is the receptacle NUMBER, 1 - indexed!)
# if no control is specified, the measurement is assumed to be stationary


def lapdReadHDF(src=None, dest=None, channel_arr = None, controls = None ):

    if controls is not None:
        motion = True
        motion_attrs = {'unit':'cm'}

        if controls[0][0] == 'NI_XZ':
            controls = None
            motion_format = 'fixed_rotation'
            with h5py.File(src, "r") as sf:
                motion_group = sf['Raw data + config/NI_XZ']
                for item in motion_group.items():
                    #Find the group, which will be the motion group
                    if isinstance(item[1], h5py.Group):
                        config_name = str(item[0])
                        #print('Motion configuration: ' +  config_name)
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

        elif controls[0][0] == 'NI_XYZ':
            motion_format = 'fixed_rotation'
            controls = None
            with h5py.File(src, "r") as sf:
                motion_group = sf['Raw data + config/NI_XYZ']
                for item in motion_group.items():
                    #Find the group, which will be the motion group
                    if isinstance(item[1], h5py.Group):
                        config_name = str(item[0])
                        #print('Motion configuration: ' +  config_name)
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


    sf = lapd.File(src, silent=True)
    

    with h5py.File(dest.file, "a") as df:
        grp = df.require_group(dest.group)
        
         #If motor drive is not included in the LAPD package, process it separately.
       
        nchan = len(channel_arr)
        for i in range(nchan):
            channel = channel_arr[i]
            data = sf.read_data(channel[2], channel[3], digitizer =channel[0],
                                adc = channel[1], add_controls=controls, 
                                silent=True)
            # Only need to do this for one channel, since we assume
            # that all channels coorespond to a single physical probe
            if i == 0:
                shp = data['signal'].shape
                nshots= shp[0]
                nti = shp[1]
                grp.require_dataset("data", (nshots, nti, nchan), np.float32, chunks=True )
                
                clock_rate = data.info['clock rate'].to(u.Hz)
                dt =  (  1.0 / clock_rate  ).to(u.s)
                time = np.arange(nti)*dt
                chan_axis = np.arange(nchan)
                shots_axis = np.arange(nshots)
                
                if controls is not None:
                    if controls[0][0] == '6K Compumotor':
                        motion_format = 'fixed_rotation'
                        pos =  data['xyz']

            
            grp['data'][:,:,i] = data['signal']
            grp['data'].attrs['unit'] = 'V'
            grp.attrs['dt'] = [s.encode('utf-8') for s in [str(dt.value), str(dt.unit)] ]
        
        
        if motion:
            grp.require_dataset('pos', (nshots, 3), np.float32)[:] = pos

            del(pos)
            for k in motion_attrs:
                grp['pos'].attrs[k] = motion_attrs[k]
                grp.attrs['motion_format'] = motion_format
            del(motion_attrs)
            
        
        dimlabels = ['shots', 'time', 'chan']
        grp['data'].attrs['dimensions'] = [s.encode('utf-8') for s in dimlabels]
        grp['data'].attrs['shape'] = np.array([nshots, nti, nchan])
        
        grp.require_dataset('shots', (nshots,), np.float32, chunks=True )[:] = shots_axis
        grp['shots'].attrs['unit'] = ''
        grp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = time.value
        grp['time'].attrs['unit'] =  str(time.unit)
        grp.require_dataset('chan', (nchan,), np.float32, chunks=True)[:] = chan_axis
        grp['chan'].attrs['unit'] = ''
        
    #Clear the LAPD HDF file from memory
    del(sf, data, time)
   
  
    


def readRunProbe( run, probe, data_dir, dest):

    csv_dir = os.path.join(data_dir, c.metadata_dir)

    run_level_attrs = csvtools.getRunLevelAttrs(csv_dir, run)
    probe_level_attrs =  csvtools.getProbeLevelAttrs(csv_dir, run, probe)
    
    attrs = {**run_level_attrs,  **probe_level_attrs}


    req_keys = ['datafile', 'digitizer', 'adc']
    missing_keys = []
    for k in req_keys:
        if not k in attrs.keys():
            missing_keys.append(k)
    if len(missing_keys) > 0:
        raise ValueError("Missing columns in csv files! The following keys were not found, and are required: " + str(missing_keys))
        
    motion_keys = ['motion_controller', 'motion_receptacle']    
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
    #Loop through the file until there are no more channels specified
    #Names must be "brd#" and "chan#"
    while True:
        brdstr = 'brd' + str(int(nchan))
        chanstr = 'chan' + str(int(nchan))
        if brdstr in attrs.keys() and chanstr in attrs.keys():
            #CHeck to make sure channel has actual non-nan values
            if not np.isnan(attrs[brdstr][0])  and not np.isnan(attrs[chanstr][0]):
                #Append the channel to the list to be extracted
                channel_arr.append( (digitizer, adc, attrs[brdstr][0], attrs[chanstr][0]) )
            nchan = nchan + 1
        else:
            break

    controls = None #Default value
    if motion:
        motion_controller = attrs['motion_controller'][0]
        motion_receptacle = attrs['motion_receptacle'][0]
        if motion_controller in ['6K Compumotor', 'NI_XZ', 'NI_XYZ']:
            controls = [(motion_controller, motion_receptacle)]


    src =  os.path.join(data_dir, c.hdf_dir, attrs['datafile'][0] +  '.hdf5')
    
    #Write the attribute dictionaries into the file
    with h5py.File(dest.file, "a") as df:
        #Delete any existing group at this location
        grp = dest.group

        #Create the new group structure.
        run_group = df.require_group('run' +str(int(run)) )
        probe_group = run_group.require_group( probe )
        dest.group = probe_group.name
        

        for k in run_level_attrs:
            if run_level_attrs[k] is not None:
                val, unit =  run_level_attrs[k]
                val = str(val).encode('utf-8')
                unit = str(unit).encode('utf-8')
                run_group.attrs.create(k,  (val, unit) ) 
                
        for k in probe_level_attrs:
            if probe_level_attrs[k] is not None:
                val, unit =  probe_level_attrs[k]
                val = str(val).encode('utf-8')
                unit = str(unit).encode('utf-8')
                probe_group.attrs.create(k,  (val, unit) ) 
     

    lapdReadHDF(src=src, dest = dest, channel_arr = channel_arr, controls = controls)
    
    

    
    
    
if __name__ == "__main__":

    data_dir = os.path.join("F:", "/LAPD_Mar2018/")

    dest = hdftools.hdfPath( r"F:/LAPD_Mar2018/RAW/test_save.hdf5")

    print('reading')
    x =  readRunProbe(102, 'PL11B', data_dir, dest)

    print('done')
    
