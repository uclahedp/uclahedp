#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program opens an LAPD HDF file (as well as metadata csvs)
and creates an ndfFile object containing the data. 

@author: peter
"""

import numpy as np
from astropy import units as u
from ndfClass import ndf
#bapsflib is available from pip. Run command 'pip papsflib' from terminal to install it
from bapsflib import lapd
from bapsflib._hdf import HDFMap




def bapsfReadHDF(filepath):
    boards =  [2,2,2]
    channels = [1,2,3]
    chan_labels = ['X', 'Y', 'Z']
    nchan = len(channels)
    
    digitizer = 'SIS crate'
    adc = 'SIS 3305'
    receptacle = 2 # This is the motion drive NUMBER, not index, so 1-indexed!!
    controls = [('6K Compumotor', receptacle)]


    f = lapd.File(filepath, silent=True)
    
    arr = []
    
    for i in range(nchan):
        data = f.read_data(boards[i], channels[i], silent=True, add_controls=controls)
        signal = data['signal'].T
        
        # Only need to do this for one channel, since we assume
        # that all channels coorespond to a single physical probe
        if i == 0:
            shotnum = data['shotnum']
            nshots = len(shotnum)
            nti = signal.shape[0]
            clock_rate = data.info['clock rate'].to(u.Hz)
            dt =  (  1.0 / clock_rate  ).to(u.s)
            
            xyz = data['xyz']
            # Get the sorted list of unique points in each direction as the axes
            # round these to .1 mm (2 decimal places) to reflect precision of LAPD drives
            xaxis = np.unique( np.around( xyz[:,0], decimals = 2 ) )*u.cm
            yaxis = np.unique( np.around( xyz[:,1], decimals = 2 ) )*u.cm
            zaxis = np.unique( np.around( xyz[:,2], decimals = 2 ) )*u.cm
            nx = len(xaxis)
            ny = len(yaxis)
            nz = len(zaxis)
            npos = nx*ny*nz
            nreps = int( nshots/npos )
        
        # Reshape the array according to the axes read in
        signal = np.reshape(signal, (nti, nx, ny, nz, nreps))
        # Temporarily store the array
        arr.append(signal)
        
    
        
    # Create an output data array with the proper form
    output = np.zeros([nti, nx, ny, nz, nreps, nchan])
    for i in range(nchan):
        output[..., i] = arr[i]
        
    # Eliminate trivial dimensions 
    output = np.squeeze(output)
    # Set units of data to volts
    output = output*u.V
    # Set datalabel
    data_label = 'Raw LAPD data'
    
    t = np.arange(nti)*dt
    rep_axis = np.arange(nreps)*u.Unit('')
    chan_axis = np.arange(nchan)*u.Unit('')
    
    axes = []
    # Assemble axis label list
    if nti > 1:
        axes.append( {'name':'t', 'axis':t} )
    if nx > 1:
        axes.append( {'name':'x', 'axis':xaxis} )
    if ny > 1:
        axes.append( {'name':'y', 'axis':yaxis} )
    if nz > 1:
        axes.append( {'name':'z', 'axis':zaxis} )
    if nreps > 1:
        axes.append( {'name':'reps', 'axis':rep_axis} )
    if nchan > 1:
        axes.append( {'name':'channels', 'axis':chan_axis} )

    obj = ndf(data=output, axes=axes, data_label = data_label)
    obj.appendLog('Created by bapsftools.bapsfReadHDF from HDF file: ' + filepath)
    
    return obj
  
    

    
    
    
if __name__ == "__main__":
    f = r"/Volumes/PVH_DATA/LAPD_Mar2018/HDF/peening056_apr05.hdf5"
    
    sfile = r"/Volumes/PVH_DATA/LAPD_Mar2018/raw/test_lapd_to_raw.h5"
    
    obj = bapsfReadHDF(f)
    obj.saveHDF( sfile)
    
