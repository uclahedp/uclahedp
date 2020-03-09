# -*- coding: utf-8 -*-
"""
process.py
@author: Peter
"""
import os
from uclahedp.load import lapd, hrr, imgdir
from uclahedp.tools import hdf, csv
from uclahedp.process import tdiode, bdot, langmuir, scope, imgseq



def process(data_dir, run, probe, 
            overwrite_raw=True, overwrite_full=True,
            csv_dir=None, hdf_dir=None, 
            img_dir=None, tdiode_hdf=None, 
            trange=None, rawsource=None):

    if csv_dir is None:
        csv_dir = os.path.join(data_dir, "METADATA")
        
    if hdf_dir is None:
        hdf_dir = os.path.join(data_dir, "HDF")
        
    if img_dir is None:
        img_dir = os.path.join(data_dir, "PIMAX")
        
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
            rawfile =  lapd.lapdToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, 
                                      trange=trange, verbose=True)
        if rawsource == 'HRR':
            print("Running hrrToRaw")
            rawfile =  hrr.hrrToRaw(run, probe[0], hdf_dir, csv_dir, rawfile, verbose=True, debug=False)
        if rawsource == 'imgdir':
            print("Running ImgDirToRaw")
            rawfile =  imgdir.imgDirToRaw(run, probe[0], img_dir, rawfile, csv_dir, verbose=True)

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
                                          calibrate=True, integrate=True, angle_correction = True,
                                          verbose=True, highfreq_calibrate=False,
                                          strict_axes=False, strict_grid=False,
                                          grid_precision = 0.25,
                                          replace_badshots = True,
                                          remove_offset=True, offset_range=(0, -50), offset_rel_t0 = [False, True])
        elif probe[1] == 'isat':
             print("Running isatRawToFull")
             fullfile = langmuir.isatRawToFull(rawfile, fullfile, tdiode_hdf=tdiode_hdf, 
                                               verbose=True, grid=True, mu=1, ti=1)
        elif probe[1] == 'vsweep':
             print("Running vsweepRawToFull")
             nfile = hdf.hdfPath(os.path.join(data_dir, "FULL", probe_string + '_ne.hdf5') )
             tfile = hdf.hdfPath(os.path.join(data_dir, "FULL", probe_string + '_te.hdf5') )
             
             fullfile = langmuir.vsweepLangmuirRawToFull(rawfile, nfile, tfile, verbose=True, grid=True, plots=True)      

        elif probe[1] == 'scope':
            print("Running scopeRawToFull")
            fullfile = scope.scopeRawToFull(rawfile, fullfile, tdiode_hdf=tdiode_hdf,
                  verbose=False, debug = False)         
        
        elif probe[1] == 'camera':
            print("Running imgseqRawToFull")
            fullfile = imgseq.imgSeqRawToFull(rawfile, fullfile)                     
                                        
        else:
            print("NO MATCHING PROBE TYPE ROUTINE EXISTS: SKIPPING!")
    else:
        print("Full file exist: skipping")
        
        

def processMany(data_dir, runs=None, probes=None, 
                overwrite_raw=True, overwrite_full=True,
                csv_dir=None, hdf_dir=None, rawsource=None,
                trange = [0,-1],
                use_tdiode='tdiode'):
    

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
        
        if (use_tdiode,'tdiode') in probelist:
            #Pop the tdiode to the front of the list, since others depend on it
            probelist.insert(0, probelist.pop(probelist.index( (use_tdiode, 'tdiode')  )))


        if len(probelist) == 0:
            print("WARNING: NO PROBES FOUND FOR RUN " + str(run))
        
        
        for probe in probelist:
            tdiode_full = hdf.hdfPath(os.path.join(data_dir, "FULL", 'run'+str(run) + '_' + use_tdiode + '_full.hdf5') )

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
                    trange=trange, rawsource=rawsource)
                


if __name__ == "__main__":
    #Windows
    #data_dir =  os.path.join("G:", "2019BIERMANN")
    #data_dir =  os.path.join("G:", "LAPD_Apr2017")
    #data_dir =  os.path.join("G:", "LAPD_Aug2016")
    #data_dir =  os.path.join("G:", "LAPD_Jan2019")
    #ata_dir =   os.path.join("G:", "LAPD_Mar2018")
    #data_dir =  os.path.join("G:", "LAPD_Jul2019")
    #data_dir =  os.path.join("G:", "LAPD_Sept2019")
    
    #OSX
    data_dir =  os.path.join("/Volumes", "PVH_DATA","2020BIERMANN")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Aug2015")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Jan2019")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Mar2018")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Jul2019")
    #data_dir =  os.path.join("/Volumes", "PVH_DATA","LAPD_Sept2019")
    
    rawsource='LAPD'
    rawsource='HRR'
    
    
    #processMany(data_dir, overwrite_raw=True, overwrite_full=True, runs=[56],probes=['tdiode'], rawsource=rawsource, trange=[0,3000])
    processMany(data_dir, overwrite_raw=True, overwrite_full=True, runs=[65],probes=['tdiode', 'PLL_B2'], use_tdiode='tdiode',  rawsource=rawsource)
    
