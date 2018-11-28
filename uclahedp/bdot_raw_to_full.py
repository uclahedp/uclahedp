#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bdot_raw_to_full.py: BDOT analysis package

Created on Wed Nov 28 13:37:21 2018

@author: peter
"""
import csvtools
import numpy as np

def bdot_raw_to_full(rawfilename, csvdir, tdiode_hdf=None):
    
    # ******
    # Load data from the raw HDF file
    # ******
    with h5py.File(rawfilename, 'r') as f:
        data = f['data']
        pos = f['pos']
        
        print(np.shape(data))
        
        run = f.attrs['run']
        probe = f.attrs['probe_name']
        
        nti = f.attrs['nti']
        npos = f.attrs['npos']
        nreps = f.attrs['nreps']
        nchan = f.attrs['nchan']
        
        daq = f.attrs['daq']
        plane = f.attrs['plane']
        drive = f.attrs['drive']
        dt = f.attrs['dt']
    
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
    atten = 10.0^(atten/20.0) # Convert from decibels
    cal = 1.0e4*(atten/gain)
    
    return


if __name__ == "__main__":
    csvdir = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/"
    rawfilename = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"
    full_file = bdot_raw_to_full(rawfilename, csvdir)
    print('Done')