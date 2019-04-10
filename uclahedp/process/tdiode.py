import h5py
import numpy as np
import os

from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util


def tdiodeRawToFull(src, dest, verbose=False, badshotratio=None,
                    fatal_badshot_percentage = None):
    
    with h5py.File(src.file, 'r') as sf:
        srcgrp= sf[src.group]
        #Get an array of all the t0 indices
        t0indarr = calcT0ind(srcgrp, verbose=verbose)
        #Get an array of all the good shots and badshots (indices)
        badshots, goodshots = findBadShots(srcgrp, verbose=verbose, 
                                           badshotratio=badshotratio,
                                           fatal_badshot_percentage=fatal_badshot_percentage)
        #Replace any bad shots with a standin avg value
        #(Later should overwrite bad shots with good neighboring shots)
        t0indarr[badshots] = int(np.median(t0indarr[goodshots]))
   
    with h5py.File(dest.file) as df:
        grp = df[dest.group]
        grp['t0indarr'] = t0indarr
        grp['badshots'] = badshots
        grp['goodshots'] = goodshots
    
    return dest



def calcT0ind(srcgrp, verbose=False):
    try:
        nshots, nti, nchan = srcgrp['data'].shape
    except KeyError:
        raise KeyError("tdiode.calcT0ind requires the data array to have an attribute 'shape'!")
        
    t0ind_array = np.zeros([nshots])
    
    tr = util.timeRemaining(nshots)
    if verbose:
        print("Calculating t0 indices")
    
    for i in range(nshots):
        #Update time remaining
        if verbose:
                tr.updateTimeRemaining(i)
        t0ind_array[i] =  np.argmax( np.gradient( srcgrp['data'][i,:,0] ) ) 

    del(nshots,nti,nchan)
    return t0ind_array.astype(int)


def calcT0(srcgrp, verbose=False):
    t0ind_array = calcT0ind(srcgrp)
    try:
        time = srcgrp['time'][:]
    except KeyError:
        raise KeyError("tdiode.calcT0 requires time array")
    return time[t0ind_array]



def findBadShots(srcgrp, verbose=False, badshotratio=None, fatal_badshot_percentage=None):
    try:
        nshots, nti, nchan = srcgrp['data'].shape
    except KeyError:
        raise KeyError("tdiode.findBadShots requires the data array to have an attribute 'shape'!")
    
    goodshots_arr = []
    badshots_arr = []
    
    if badshotratio is None:
        badshotratio = 10
        
    if fatal_badshot_percentage is None:
        fatal_badshot_percentage= 0.2
    
    tr = util.timeRemaining(nshots)
    if verbose:
        print("Identifying bad shots")
        
    for i in range(nshots):
        #Update time remaining
        if verbose:
                tr.updateTimeRemaining(i)
        #TODO: trying using the mean of the last 500 points rather than the median as the reference
        max_median_ratio = np.max( srcgrp['data'][i,:,0]) / np.mean(srcgrp['data'][i,-500:,0])
        #This defines a 'bad shot' where the laser diode was indistinct,
        #indicating a possible misfire
        if (max_median_ratio < badshotratio):
            badshots_arr.append(i)
        else:
            goodshots_arr.append(i)
            
    print("Found " + str(len(badshots_arr)) + ' bad shots')
    
    if len(badshots_arr)/nshots > fatal_badshot_percentage:
        raise ValueError("Lots of bad shots found! Bad sign! Aborting.")
            
                
    return badshots_arr, goodshots_arr




    



if __name__ == "__main__":
    #src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode_raw.hdf5"))
    #dest = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode_full.hdf5"))
    
    probe = 'tdiode'
    run = 6
    
    src = hdftools.hdfPath( '/Volumes/PVH_DATA/2019BIERMANN/RAW/' + 'run' + str(run) + '_' + probe + '_raw.hdf5')
    dest = hdftools.hdfPath('/Volumes/PVH_DATA/2019BIERMANN/FULL/' + 'run' + str(run) + '_' + probe + '_full.hdf5')


    #Delete the output file if it already exists
    try:
        os.remove(dest.file)
    except FileNotFoundError:
        pass
    
    #print(calcT0(src)[0:20] )
    #print(findBadShots(src) )
    tdiodeRawToFull(src, dest, verbose=True)