#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bdot.py: BDOT analysis package

Created on Wed Nov 28 13:37:21 2018

@author: peter
"""

# TODO 
# Program should check for four numbers in the spreadsheet: x, y, z, roll
# numbers from the spreadsheet should take priority. Otherwise, use values from hdf5
# Also allow a keyword argument to the function for each of these to allow setting manually


import csvtools
import numpy as np
import os
from shutil import copyfile
import h5py

def bdot_raw_to_full(rawfilepath, csvdir, tdiode_hdf=None, fullfilepath=None):
    """ Integrates bdot data, calibrates output using information about the probe.
        Corrects for probe angle based on which drive is being used.

    Parameters
    ----------
        rawfilepath: str
            Path string to a raw hdf5 file containing bdot data
            
        csvdir: str
            Path to the directory that contains the csv metadata files
            The names of the files in that directory are hardcoded in here...
            
        tdiode_hdf: str
            Path to a tdiode full hdf5 file containing diode t0 times and
            an array of boolean 'bad shot' flags.
            
        fullfilepath: str
            Path at which to store the full hdf file once it is created

    Returns
    -------
        fullfilepath : str
            Path string for the newly created full hdf5 file
    """ 
    # ******
    # Load data from the raw HDF file
    # ******
    with h5py.File(rawfilepath, 'r') as f:
        data = f['data'][:]
        pos = f['pos'][:]
    
        run = f.attrs['run']
        probe = f.attrs['probe_name']
        
        nti = f.attrs['nti']
        npos = f.attrs['npos']
        nreps = f.attrs['nreps']
        nchan = f.attrs['nchan']
        
        drive = f.attrs['drive']
        dt = f.attrs['dt'] # dt MUST be in s for this algorithm to work...

    
    # ******
    # Load data from bdot runs csv
    # ******
    bdot_runs_csv = csvdir + "bdot_runs.csv"
    bdot_runs_csv = csvtools.opencsv(bdot_runs_csv)
    
    xatten = csvtools.findvalue(bdot_runs_csv, 'xatten', run=run, probe=probe)
    yatten = csvtools.findvalue(bdot_runs_csv, 'yatten', run=run, probe=probe)
    zatten = csvtools.findvalue(bdot_runs_csv, 'zatten', run=run, probe=probe)
    atten = np.array([xatten, yatten, zatten])
    
    xpol = csvtools.findvalue(bdot_runs_csv, 'xpol', run=run, probe=probe)
    ypol = csvtools.findvalue(bdot_runs_csv, 'ypol', run=run, probe=probe)
    zpol = csvtools.findvalue(bdot_runs_csv, 'zpol', run=run, probe=probe)
    pol = np.array([xpol, ypol, zpol])
    
    probe_origin_x = csvtools.findvalue(bdot_runs_csv, 'probe_origin_x', run=run, probe=probe)
    probe_origin_y = csvtools.findvalue(bdot_runs_csv, 'probe_origin_y', run=run, probe=probe)
    probe_origin_z = csvtools.findvalue(bdot_runs_csv, 'probe_origin_z', run=run, probe=probe)
    probe_origin = np.array([probe_origin_x, probe_origin_y, probe_origin_z])

    gain = csvtools.findvalue(bdot_runs_csv, 'gain', run=run, probe=probe)
    probe_rot = csvtools.findvalue(bdot_runs_csv, 'probe_rot', run=run, probe=probe)

    # ******
    # Load data from bdot probes csv
    #*******
    bdot_probes_csv = csvdir + "bdot_probes.csv"
    bdot_probes_csv = csvtools.opencsv(bdot_probes_csv)
    
    xarea = csvtools.findvalue(bdot_probes_csv, 'xarea',  probe=probe)
    yarea = csvtools.findvalue(bdot_probes_csv, 'yarea',  probe=probe)
    zarea = csvtools.findvalue(bdot_probes_csv, 'zarea',  probe=probe)
    area = np.array([xarea, yarea, zarea])
    nturns = csvtools.findvalue(bdot_probes_csv, 'num_turns',  probe=probe)
    
    # ******
    # Load data from main runs csv
    # ******
    main_runs_csv = csvdir + "main_runs.csv"
    main_runs_csv = csvtools.opencsv(main_runs_csv)
    
    target_xpos = csvtools.findvalue(main_runs_csv, 'target_xpos', run=run)
    target_ypos = csvtools.findvalue(main_runs_csv, 'target_ypos', run=run)
    target_zpos = csvtools.findvalue(main_runs_csv, 'target_zpos', run=run)
    target_pos = np.array([target_xpos, target_ypos, target_zpos])
    
    # TODO: add this functionality
    if tdiode_hdf is not None:
        print("Handle tdiode correction here...")
        
    # Create an array of calibration coefficents
    atten = np.power([10,10,10], atten/20.0) # Convert from decibels
    # Required input units
    # dt -> s
    # dt is already in s
    # area : mm^2 -> m^2
    area = area*1e-6
    cal = 1.0e4*dt*atten/gain/(nturns*area)


    #Integrate the data
    bx = np.cumsum(data[:, :, :, 0], axis=0)
    by = np.cumsum(data[:, :, :, 1], axis=0)
    bz = np.cumsum(data[:, :, :, 2], axis=0)
    
    bx = cal[0]*pol[0]*bx
    by = cal[1]*pol[1]*by
    bz = cal[2]*pol[2]*bz
    
    
    # Correct for probe rotation (generally accidental...)
    # This is rotation about the probe's main (x) axis
    if probe_rot is not 0:
        probe_rot = np.deg2rad(probe_rot)
        by = by*np.cos(probe_rot) - bz*np.sin(probe_rot)
        bz = bz*np.cos(probe_rot) + by*np.sin(probe_rot)

    #Reassemble the data array, correcting for angles due to the drive
    
    if drive in ['none', 'cartesian_xyz' ]:
        data[:, :, :, 0] = bx
        data[:, :, :, 1] = by
        data[:, :, :, 2] = bz
    elif drive in ['xy','polar_xy']:
        # Calculate the angle made by the probe shaft at each pos
        angle = np.arctan(pos[1,:]/pos[0,:])
        # The remainder of this mess creates a matrix ready for multiplication
        theta = np.outer(np.ones(nti), angle)
        theta = theta.flatten()
        theta = np.outer(theta,np.ones(nreps))
        theta = np.reshape( theta.flatten(), [nti,npos, nreps])
    
        data[:, :, :, 0] = bx*np.cos(theta) - by*np.sin(theta)
        data[:, :, :, 1] = by*np.cos(theta) + bx*np.sin(theta)
        data[:, :, :, 2] = bz
    elif drive in ['xz','polar_xz']:
        # Calculate the angle made by the probe shaft at each pos
        angle = np.arctan(pos[2, :]/pos[0, :])
        # The remainder of this mess creates a matrix ready for multiplication
        theta = np.outer(np.ones(nti), angle)
        theta = theta.flatten()
        theta = np.outer(theta,np.ones(nreps))
        theta = np.reshape( theta.flatten(), [nti,npos, nreps])
        
        data[:, :, :, 0] = bx*np.cos(theta) - bz*np.sin(theta)
        data[:, :, :, 1] = by
        data[:, :, :, 2] = bz*np.cos(theta) + bx*np.sin(theta)
    else:
        print("Invalid drive type: " + drive)
        return None
    
    
    # Shift the position array to be relative to TCC
    pos[0,:] = pos[0,:] - target_pos[0]
    pos[1,:] = pos[1,:] - target_pos[1]
    pos[2,:] = pos[2,:] - target_pos[2]
    


    #If necessary, come up with a new filename
    if fullfilepath is None:
        rawdir = os.path.dirname(rawfilepath) + '/'
        fullfilepath = rawdir + 'run' + str(run) + '_' + str(probe) + '_full.h5'
        
    copyfile(rawfilepath, fullfilepath)
    
    #Change the relevant variables in the previous HDF file
    with h5py.File(fullfilepath, 'a') as f2:
        f2['data'][...] = data # The [...] is essential for overwriting data, but I don't understand what it does...
        f2["data"].attrs['unit'] = 'G'
        f2['pos'][...] = pos
        f2.attrs['data_type'] = 'full'
        f2.attrs['chan_labels'] =  [s.encode('utf-8') for s in ['BX', 'BY', 'BZ']] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289 

    return fullfilepath


if __name__ == "__main__":
    csvdir = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/"
    rawfilepath = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    #rawfilename = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_pos_raw.h5"
    full_filepath = bdot_raw_to_full(rawfilepath, csvdir)
    print('Done')