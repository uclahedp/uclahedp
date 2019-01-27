import h5py
import numpy as np

import os
import hdftools

#TODO: Write a tdiode raw_to_full program so other probes don't have to regen
#this info for each probe



def calcT0ind(src):
    with h5py.File(src.file, 'r') as sf:
        grp = sf[src.group]

        try:
            nshots, nti, nchan = grp['data'].attrs['shape']
        except KeyError:
            raise KeyError("tdiode.calcT0ind requires the data array to have an attribute 'shape'!")
            
        t0ind_array = np.zeros([nshots])
        for i in range(nshots):
            t0ind_array[i] =  np.argmax( np.gradient( grp['data'][i,:,0] ) ) 

    del(nshots,nti,nchan)
    return t0ind_array.astype(int)


def calcT0(src):
    t0ind_array = calcT0ind(src)
    with h5py.File(src.file, 'r') as sf:
        grp = sf[src.group]
        try:
            time = grp['time'][:]
        except KeyError:
            raise KeyError("tdiode.calcT0 requires time array")
    return time[t0ind_array]



def findBadShots(src):
    with h5py.File(src.file, 'r') as sf:
        grp = sf[src.group]
        try:
            nshots, nti, nchan = grp['data'].attrs['shape']
        except KeyError:
            raise KeyError("tdiode.findBadShots requires the data array to have an attribute 'shape'!")
        
        goodshots_arr = []
        badshots_arr = []
        for i in range(nshots):
            max_median_ratio = np.max( grp['data'][i,:,0]) / np.median(grp['data'][i,:,0])
            #This defines a 'bad shot' where the laser diode was indistinct,
            #indicating a possible misfire
            if max_median_ratio < 10:
                badshots_arr.append(i)
            else:
                goodshots_arr.append(i)
                
        return badshots_arr, goodshots_arr




    



if __name__ == "__main__":
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_save.hdf5"), 'run102/tdiode')
    
    print(calcT0(src)[0:20] )
    print(findBadShots(src) )