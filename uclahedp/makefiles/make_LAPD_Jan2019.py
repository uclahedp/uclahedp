# -*- coding: utf-8 -*-
"""
Runscript for LAPD_Jan2019 experiment
@author: Peter
"""
import os
from uclahedp import lapdtools, hdftools, csvtools, tdiode, bdot


#Set this flag to delete files that already exist and process them again
#Otherwise, the program will assume they are already ready
overwrite = True

#Windows
#data_dir =  os.path.join("F:", "LAPD_Jan2019")
#OSX
data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Jan2019")

#Define some other paths
csv_dir = os.path.join(data_dir, "METADATA")
hdf_dir = os.path.join(data_dir, "HDF")

#Generate a list of all available runs OR set your own runlist
runlist = csvtools.getRunList(csv_dir)
runlist = [18]


#Set a list of probes to run here
#If its empty, ALL OF THE PROBES WILL BE RUN
probes = ['tdiode']


for run in runlist:
    
    print("Processing run: " + str(run))
    
    probelist = csvtools.getProbeList(csv_dir, run)
    #If a tdiode exists, bring it to the front of the list so it gets called first
    try:
        probelist.insert(0, probelist.pop(probelist.index( ('tdiode', 'tdiode')  )))
    except ValueError:
        pass
    
    if len(probelist) == 0:
        print("WARNING: NO PROBES FOUND FOR RUN " + str(run))
    
    
    for probe in probelist:
        
        #If the probes array is set but doesn't include this probe,
        #skip it
        if len(probes) != 0 and probe[0] not in probes:
            print("SKIPPING: " + probe[0] + ', NOT IN PROBE LIST!')
            continue
        

        print("Processing probe: " + str(probe))
        probe_string = 'run'+str(run) + '_' + probe[0]
        
        
        if probe[1] == 'tdiode' or probe[1] == 'unknown':
            full_dir = os.path.join(data_dir, "TDIODE")
        elif probe[1] == 'bdot':
            full_dir = os.path.join(data_dir, "BDOT")
            
            
        rawfile = hdftools.hdfPath(os.path.join(data_dir, "RAW", probe_string + '_raw.hdf5') )
        fullfile = hdftools.hdfPath(os.path.join(full_dir, probe_string + '_full.hdf5') )
        tdiode_full = hdftools.hdfPath(os.path.join(data_dir, "TDIODE", 'run'+str(run) + '_tdiode_full.hdf5') )
        
        
        #Make the raw file
        if os.path.exists(rawfile.file):
            if overwrite:
                os.remove(rawfile.file)
        if not os.path.exists(rawfile.file):
            print("Running lapdToRaw")
            rawfile =  lapdtools.lapdToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, verbose=True)
        else:
            print("Raw file exists: skipping")
            
            
            
        #Make the full file
        if os.path.exists(fullfile.file):
            if overwrite:
                os.remove(fullfile.file)
        if not os.path.exists(fullfile.file):
            if probe[1] == 'tdiode':
                print("Running tdiodeRawToFull")
                fullfile = tdiode.tdiodeRawToFull(rawfile, fullfile, verbose=True)
            elif probe[1] == 'bdot':
                print("Running bdotRawToFull")
                fullfile = bdot.bdotRawToFull(rawfile, fullfile, tdiode_hdf=tdiode_full, grid=True, verbose=True)
            else:
                print("NO MATCHING PROBE TYPE ROUTINE EXISTS: SKIPPING!")
        else:
            print("Full file exists: skipping")
                

