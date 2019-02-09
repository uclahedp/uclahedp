import h5py
import numpy as np
import os

from uclahedp import hdftools, util

#TODO: Write a tdiode raw_to_full program so other probes don't have to regen
#this info for each probe


def tdiodeRawToFull(src, dest, verbose=False):
    with h5py.File(src.file, 'r') as sf:
        srcgrp= sf[src.group]
        #Get an array of all the t0 indices
        t0indarr = calcT0ind(srcgrp, verbose=verbose)
        #Get an array of all the good shots and badshots (indices)
        badshots, goodshots = findBadShots(srcgrp, verbose=verbose)
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



def findBadShots(srcgrp, verbose=False):
    try:
        nshots, nti, nchan = srcgrp['data'].shape
    except KeyError:
        raise KeyError("tdiode.findBadShots requires the data array to have an attribute 'shape'!")
    
    goodshots_arr = []
    badshots_arr = []
    
    tr = util.timeRemaining(nshots)
    if verbose:
        print("Identifying bad shots")
        
    for i in range(nshots):
        #Update time remaining
        if verbose:
                tr.updateTimeRemaining(i)
        max_median_ratio = np.max( srcgrp['data'][i,:,0]) / np.median(srcgrp['data'][i,:,0])
        #This defines a 'bad shot' where the laser diode was indistinct,
        #indicating a possible misfire
        if max_median_ratio < 10:
            badshots_arr.append(i)
        else:
            goodshots_arr.append(i)
                
    return badshots_arr, goodshots_arr




    



if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode_raw.hdf5"))
    dest = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_tdiode_full.hdf5"))
    #print(calcT0(src)[0:20] )
    #print(findBadShots(src) )
    tdiodeRawToFull(src, dest, verbose=True)