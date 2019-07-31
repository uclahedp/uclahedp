# -*- coding: utf-8 -*-
"""
process.py
@author: Peter
"""
import os
from uclahedp.load import lapd, hrr
from uclahedp.tools import hdf, csv
from uclahedp.process import tdiode, bdot, langmuir, monochromator



def process(data_dir, run, probe, 
            overwrite_raw=True, overwrite_full=True,
            csv_dir=None, hdf_dir=None, tdiode_hdf=None, rawsource=None):

    if csv_dir is None:
        csv_dir = os.path.join(data_dir, "METADATA")
        
    if hdf_dir is None:
        hdf_dir = os.path.join(data_dir, "HDF")
        
    if rawsource is None:
        rawsource = 'LAPD'
        
    probe_string = 'run'+str(run) + '_' + probe[0]

    rawfile = hdf.hdfPath(os.path.join(data_dir, "RAW", probe_string + '_raw.hdf5') )
    fullfile = hdf.hdfPath(os.path.join(data_dir, "FULL", probe_string + '_full.hdf5') )
    
    #Make the directory if necessary
    os.makedirs(os.path.dirname(rawfile.file), exist_ok=True)
    os.makedirs(os.path.dirname(fullfile.file), exist_ok=True)

    #Make the raw file
    if os.path.exists(rawfile.file):
        if overwrite_raw:
            os.remove(rawfile.file)
    if not os.path.exists(rawfile.file):
        if rawsource == 'LAPD':
            print("Running lapdToRaw")
            rawfile =  lapd.lapdToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, verbose=True)
        if rawsource == 'HRR':
            print("Running lapdToRaw")
            rawfile =  hrr.hrrToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, verbose=True)
    else:
        print("Raw file exist: skipping")
    
    print(tdiode_hdf)
    #Make the full file
    if os.path.exists(fullfile.file):
        if overwrite_full:
            os.remove(fullfile.file)
    if not os.path.exists(fullfile.file):
        if probe[1] == 'tdiode':
            print("Running tdiodeRawToFull")
            fullfile = tdiode.tdiodeRawToFull(rawfile, fullfile, verbose=True, fatal_badshot_percentage=.2, badshotratio=10)
        elif probe[1] == 'bdot':
            print("Running bdotRawToFull")
            fullfile = bdot.bdotRawToFull(rawfile, fullfile, tdiode_hdf=tdiode_hdf, grid=True, 
                                          verbose=True, highfreq_calibrate=True,
                                          remove_offset=True, offset_range=(0, -1), offset_rel_t0 = [False, True])
        elif probe[1] == 'isat':
             print("Running isatRawToFull")
             fullfile = langmuir.isatRawToFull(rawfile, fullfile, tdiode_hdf=tdiode_hdf, 
                                               verbose=True, grid=True, mu=4, ti=1)
        elif probe[1] == 'vsweep':
             print("Running vsweepRawToFull")
             nfile = hdf.hdfPath(os.path.join(data_dir, "FULL", probe_string + '_density.hdf5') )
             tfile = hdf.hdfPath(os.path.join(data_dir, "FULL", probe_string + '_temperature.hdf5') )
             fullfile = langmuir.vsweepLangmuirRawToFull(rawfile, nfile, tfile, verbose=True, grid=True)      

        elif probe[1] == 'monochromator':
            print("Running monochromatorRawToFull")
            fullfile = monochromator.monochromatorRawToFull(rawfile, fullfile, port=14, tdiode_hdf=tdiode_hdf,
                  verbose=False, debug = False)                           
                                        
        else:
            print("NO MATCHING PROBE TYPE ROUTINE EXISTS: SKIPPING!")
    else:
        print("Full file exist: skipping")
        
        

def processMany(data_dir, runs=None, probes=None, 
                overwrite_raw=True, overwrite_full=True,
                csv_dir=None, hdf_dir=None, rawsource=None):
    

    if csv_dir is None:
        csv_dir = os.path.join(data_dir, "METADATA")
        
    if hdf_dir is None:
        hdf_dir = os.path.join(data_dir, "HDF")
        
    if runs is None:
        runlist = csv.getRunList(csv_dir)
    else:
        runlist = runs
        


    for run in runlist:
        
        print("Processing run: " + str(run))
        
        #Load the probe list
        probelist = csv.getProbeList(csv_dir, run)
        #If a tdiode exists, bring it to the front of the list so it gets called first
        
        if ('tdiode','tdiode') in probelist:
            #Pop the tdiode to the front of the list, since others depend on it
            probelist.insert(0, probelist.pop(probelist.index( ('tdiode', 'tdiode')  )))


        if len(probelist) == 0:
            print("WARNING: NO PROBES FOUND FOR RUN " + str(run))
        
        
        for probe in probelist:
            
            if probe[0] ==  'monochromator':
                tdiode_full = hdf.hdfPath(os.path.join(data_dir, "FULL", 'run'+str(run) + '_scope_tdiode_full.hdf5') )
            else:
                tdiode_full = hdf.hdfPath(os.path.join(data_dir, "FULL", 'run'+str(run) + '_tdiode_full.hdf5') )
                
                
            print(tdiode_full.file)  
            file_exists = os.path.exists(tdiode_full.file)
            if file_exists:
                print("Using tdiode file: " + str(tdiode_full.file))
            else:
                print("NO TDIODE FOUND: WILL NOT CORRECT FOR TIMING!")
                tdiode_full = None
                
                
            

            #If the probes array is set but doesn't include this probe,
            #skip it
            if probes is not None and probe[0] not in probes:
                print("SKIPPING: " + probe[0] + ', NOT IN PROBE LIST!')
                continue
            print("Processing probe: " + str(probe))
            process(data_dir, run, probe, 
                    overwrite_raw=overwrite_raw, overwrite_full=overwrite_full,
                    csv_dir=csv_dir, hdf_dir=hdf_dir, tdiode_hdf=tdiode_full,
                    rawsource=rawsource)
                


if __name__ == "__main__":
    #Windows
    #data_dir =  os.path.join("F:", "2019BIERMANN")
    #data_dir =  os.path.join("F:", "LAPD_Apr2017")
    #data_dir =  os.path.join("F:", "LAPD_Jan2019")
    #data_dir =  os.path.join("F:", "LAPD_Mar2018")
    data_dir =  os.path.join("G:", "LAPD_Jul2019")
    
    #OSX
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","2019BIERMANN")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Aug2015")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Jan2019")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Mar2018")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Jul2019")
    
    rawsource='LAPD'
    #rawsource='HRR'
    
    
    #'tdiode', 'LAPD3', 'LAPD_C2'
    #10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
    
    processMany(data_dir, overwrite_raw=True, overwrite_full=True, runs=[28,29,30,31,32,33],probes=['tdiode', 'LAPD3', 'LAPD_C2'], rawsource=rawsource) 
    #processMany(data_dir, overwrite=False, runs=[18], probes=['LAPD_C6'], rawsource=rawsource) 
    
    
