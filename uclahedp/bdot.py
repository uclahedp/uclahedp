#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bdot.py: BDOT analysis package

Created on Wed Nov 28 13:37:21 2018

@author: peter
"""

import csvtools
import numpy as np
import os
from shutil import copyfile
import h5py
import hdftools

def bdot_raw_to_full(src, dest, tdiode_hdf=None):
    """ Integrates bdot data, calibrates output using information about the probe.
        Corrects for probe angle based on which drive is being used.

    Parameters
    ----------
        src: hdfPath object
            Path string to a raw hdf5 file containing bdot data
            
        dest: hdfPath object
            Path string to location processed bdot data should be written out

        tdiode:  hdfPath object
            Path to a tdiode full hdf5 dataset containing diode t0 times and
            an array of boolean 'bad shot' flags.


    Returns
    -------
       None
    """ 
    # ******
    # Load data from the raw HDF file
    # ******
    print(src.file)
    with h5py.File(src.file, 'r') as sf:
        
        kdict = {}
        srcgrp = sf[src.group]
        
        
        #Check for keys always required by this function
        req_keys = ['brd1','brd2','brd3', 'chan1','chan2', 'chan3', 
                    'xarea', 'yarea', 'zarea',
                    'xatten', 'yatten', 'zatten', 'gain',
                    'xpol', 'ypol', 'zpol', 'roll', 
                    'probe_origin_x', 'probe_origin_y', 'probe_origin_z',
                    'dt']
       
        motion = None
        if  'pos' in srcgrp:
            req_keys = req_keys + ['motion_format']
            #If pos array exists, there are keywords required for that too.
            if srcgrp.attrs['motion_format'] == 'fixed_rotation':
                req_keys = req_keys + ['rot_center_x', 'rot_center_y', 'rot_center_z']
            else:
                raise ValueError("Motion format unrecognized: " + str(srcgrp['pos'].attrs['motion_format']) )
            motion = srcgrp.attrs['motion_format']
        else:
            #If no position information is given, a single explicit position
            #is required. 
            req_keys = req_keys + ['xpos', 'ypos', 'zpos']
            
            
        print(req_keys)
        #Process the required keys, throwing an error if any cannot be found
        missing_keys = []
        for k in req_keys:
            if k in srcgrp.attrs:
                kdict[k] = ( csvtools.fixType(srcgrp.attrs[k][0]), srcgrp.attrs[k][1])
            else:
                missing_keys.append(k)
        if len(missing_keys) > 0:
            raise ValueError("Missing required keys! ->" + str(missing_keys))
        
        
       
        

                
                
        nshots, nti, nchan = srcgrp['data'].shape
        
        with h5py.File(dest.file, 'a') as df:
            #Clear group if it already exists
            try:
                del(df[dest.group])
            except KeyError:
                pass
            
            destgrp = df.require_group(dest.group)
            destgrp.create_dataset('data', (nshots, nti, nchan))
            
            
            atten = np.array([kdict['xatten'][0],kdict['yatten'][0],kdict['zatten'][0]])
            # Create an array of calibration coefficents
            if kdict['xatten'][1] == 'dB':
                print("Converting dB to x")
                atten = np.power([10,10,10], atten/20.0) # Convert from decibels
            
            
            xarea = (kdict['xarea'][0])*u.Unit(kdict['xarea'][1]).to(u.m ** 2)
            print(xarea)
            #area = np.array([kdict['yarea'][0],kdict['zarea'][0]])
            
            
            # Required input units
            # dt -> s
            # dt is already in s
            # area : mm^2 -> m^2
            area = area*1e-6
            cal = 1.0e4*dt*atten/gain/(nturns*area)
            
            #Chunking data processing loop limits memory usage
            for i in range(nshots):
                
                bx = srcgrp['data'][i,:, 0]
                by = srcgrp['data'][i,:, 1]
                bz = srcgrp['data'][i,:, 2]
                
                if motion == 'fixed_rotation':
                    x,y,z = srcgrp['pos'][i, :]
                    rx, ry, rz = kdict['rot_center_x'][0],kdict['rot_center_y'][0],kdict['rot_center_z'][0]
                    pitch = np.arctan( (y-ry) / (x-rx) ) 
                    yaw = np.arctan( (z-rz) / (x-rx) ) 
                    
                    roll, unit = kdict['roll']
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
    

            
            del(bx,by,bz)
            
            if motion == 'fixed_rotation':
                del(x,y,z,rx,ry,rz,roll, pitch, yaw)
                
        


       
        pos = f['pos'][:]
    
        run = f.attrs['run']
        probe = f.attrs['probe_name']
        
        nti = f.attrs['nti']
        npos = f.attrs['npos']
        nreps = f.attrs['nreps']
        nchan = f.attrs['nchan']
        
        drive = f.attrs['drive']
        dt = f.attrs['dt'] # dt MUST be in s for this algorithm to work...

        
    


    #Integrate the data
    bx = np.cumsum(data[:, :, :, 0], axis=0)
    by = np.cumsum(data[:, :, :, 1], axis=0)
    bz = np.cumsum(data[:, :, :, 2], axis=0)
    
    bx = cal[0]*pol[0]*bx
    by = cal[1]*pol[1]*by
    bz = cal[2]*pol[2]*bz
    
    



    


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
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_save.hdf5"), 'run102/PL11B')
    dest = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_save_full.hdf5"), 'run102/PL11B')
    
    #rawfilename = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_pos_raw.h5"
    full_filepath = bdot_raw_to_full(src, dest)
    print('Done')
