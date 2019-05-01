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
        
        t = srcgrp['time'][:]
        nshots = srcgrp['shots'].shape[0]
        nti = srcgrp['time'].shape[0]
        nchan = srcgrp['chan'].shape[0]
        
        #Create the destination file directory if necessary
        hdftools.requireDirs(dest.file)
        
        with h5py.File(dest.file) as df:
            destgrp = df[dest.group]
            destgrp['t0indarr'] = t0indarr
            destgrp['badshots'] = badshots
            destgrp['goodshots'] = goodshots
            
    
            #Apply the tdiode correction to the data, as a check
            #If this is working correctly, the full tdiode files will show
            #the tdiode's all lined up...
            min_t0ind = np.min(t0indarr[goodshots])
            max_t0shift = np.max(t0indarr[goodshots]) - min_t0ind
            #Compute new nti
            nti = nti - max_t0shift
            
            destgrp.require_dataset('data', (nshots, nti, nchan), np.float32, chunks=(1, np.min([nti, 20000]), 1), compression='gzip')
            destgrp['data'].attrs['dimensions'] = srcgrp['data'].attrs['dimensions']
            destgrp['data'].attrs['unit'] = srcgrp['data'].attrs['unit']
            
            for i in range(0, nshots):
                ta = t0indarr[i] - min_t0ind
                tb = ta + nti
                destgrp['data'][i, :, :] = srcgrp['data'][i, ta:tb, :]
            

            destgrp.require_dataset('time', (nti,), np.float32, chunks=True)[:] = t[0:nti] - t[min_t0ind]
            destgrp['time'].attrs['unit'] = srcgrp['time'].attrs['unit']
            
            destgrp.require_dataset('chan', (nchan,), np.int32, chunks=True)[:] = srcgrp['chan'][:]
            destgrp['chan'].attrs['unit'] = srcgrp['chan'].attrs['unit']
            destgrp.require_dataset('shots', (nshots,), np.int32, chunks=True)[:] = srcgrp['shots'][:]
            destgrp['shots'].attrs['unit'] = srcgrp['shots'].attrs['unit']
            
    
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