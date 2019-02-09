# -*- coding: utf-8 -*-
"""
Runscript for LAPD_Jan2019 experiment
@author: Peter
"""
import os
from uclahedp import lapdtools, hdftools, csvtools, tdiode, bdot


def fullDir(probe):
    if probe[1] == 'tdiode':
        full_dir = os.path.join(data_dir, "TDIODE")
    elif probe[1] == 'bdot':
        full_dir = os.path.join(data_dir, "BDOT")  
    elif probe[1] == 'langmuir':
        full_dir = os.path.join(data_dir, "LANGMUIR")
    else:
        print("PROBE NOT FOUND: " + probe)
    return full_dir



overwrite = True

data_dir =  os.path.join("F:", "LAPD_Jan2019")
csv_dir = os.path.join(data_dir, "METADATA")
hdf_dir = os.path.join(data_dir, "HDF")




run = 29
probe = ('LAPD7', 'bdot')
#probe = ('tdiode', 'tdiode')



    

runlist = csvtools.getRunList(csv_dir)



probelist = csvtools.getProbeList(csv_dir, run)
#If a tdiode exists, bring it to the front of the list so it gets called first
try:
    probelist.insert(0, probelist.pop(probelist.index( ('tdiode', 'tdiode')  )))
except ValueError:
    pass


full_dir = fullDir(probe)
probe_string = 'run'+str(run) + '_' + probe[0]

rawfile = hdftools.hdfPath(os.path.join(data_dir, "RAW", probe_string + '_raw.hdf5') )
tdiode_full = hdftools.hdfPath(os.path.join(data_dir, "TDIODE", 'run'+str(run) + '_tdiode_full.hdf5') )
fullfile = hdftools.hdfPath(os.path.join(full_dir, probe_string + '_full.hdf5') )


#Make the raw file
if os.path.exists(rawfile.file):
    if overwrite:
        os.remove(rawfile.file)
if not os.path.exists(rawfile.file):
    print("Running lapdToRaw")
    rawfile =  lapdtools.lapdToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, verbose=True)
else:
    print("Raw file exist: skipping")
    
    
    
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
    print("Full file exist: skipping")
        

