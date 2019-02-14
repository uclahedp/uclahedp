#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 09:20:38 2019

@author: peter
"""
import h5py
import numpy as np

from uclahedp import hdftools, util

def validDataset(grp):
    """
    Determines whether a given HDF5 dataset group conforms to the HEDP group's 
    standard (as explained in comments below)
    """
  
    #The dataset group must be a valid HDF5 group
    if not isinstance(grp, h5py.Group ):
        raise hdfDatasetFormatError("Group is not a valid HDF5 Group!"
                                    + str(grp))
    
    #Group must have a dataset named 'data'
    if 'data' not in grp:
        raise hdfDatasetFormatError("Group must have member dataset 'data'!")
    
    #These attributes must exist on the 'data' dataset
    req_keys = ['dimensions', 'unit']
    for k in req_keys:
        if k not in grp['data'].attrs.keys():
            raise hdfDatasetFormatError("'data' must have attribute: " 
                                    + str(k) + "!")
            
    dims = grp['data'].attrs['dimensions']
    shape = grp['data'].shape
    
    for i in range(len(shape)):
        #Each dimension of the dataset must have an axes dataset
        if not dims[i] in grp:
            raise hdfDatasetFormatError("Group must have an axis for each " +
                                        "dimension, and is missing: " + 
                                        str(dims[i]) + "!")
        axis = grp[dims[i]]
        
        #Axes must have the same length as the associated dataset dimension
        if not len(axis) == shape[i]:
            raise hdfDatasetFormatError("Axis length must match associated "+
                                        "dimension, and this one does not: " +
                                        str(axis) + "!")
        #Axes must have 'unit' attribute
        if not 'unit' in axis.attrs.keys():
            raise hdfDatasetFormatError("Axis must have a 'unit' attribute, " +
                                        "and this one does not: " +
                                        str(axis) + "!")

    return True


class hdfDatasetFormatError(Exception):
    """
    This HDF5 dataset does not conform to the HEDP group's standard data layout
    """
    pass



def avgDim(src, dest, ax, delsrc=False, verbose=False):
    chunked_array_op(src, dest, ax, 'avgDim', delsrc=delsrc, verbose=verbose)
    
def avgDimOp(src_dset, dest_dset, sl, axind, args):
    s = src_dset[tuple(sl)]
    s = np.mean(s, axis=axind, keepdims=True)
    sl[axind] = slice(None, None, None)  # This dim is 1 in the dest dataset
    dest_dset[tuple(sl)] = s
    
    
def thinPick(src, dest, ax, step=None, delsrc=False, verbose=False):
    if step is None:
        step = 10
    else:
        step = int(step)    
    chunked_array_op(src, dest, ax, 'thinPick', delsrc=delsrc, verbose=verbose, step=step)

def thinPickOp(src_dset, dest_dset, sl, axind, args):
    #Thin by plucking out entries
    sl[axind] = slice(None, None, args['step'])
    s = src_dset[tuple(sl)]
    sl[axind] = slice(None, None, None)
    print(s.shape)
    dest_dset[tuple(sl)] = s


        



def chunked_array_op(src, dest, ax, oplabel, delsrc=False, verbose=False, **args):
    
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        
        validDataset(srcgrp)
        
        oldshape = list( srcgrp['data'].shape )
        ndim = len(oldshape)
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        
        try:
            axind = dimlabels.index(ax)
        except ValueError as e:
            print("ERROR: Couldn't find ax in dimensions array!: " + str(ax) )
            raise(e)
        
        #Decide on a chunking axis
        #Get a list of the axes indices ordered by chunk size, largest to smallest
        chunks =  np.flip( np.argsort( srcgrp['data'].chunks ) )
        
        #Chose the largest one that ISN'T the chosen axis
        chunkax = chunks[0]
        if chunkax == axind:
            chunkax = chunks[1]
            
        chunksize = 1000
        nchunks = int( np.ceil(oldshape[chunkax] / chunksize))

        #Make alterations to the shape for the new dataset
        #This will obviously depend on what operation is being performed
        newshape = np.copy(oldshape)
        if oplabel == 'avgDim':
            newshape[axind] = 1
        elif oplabel == 'thinPick':
            newshape[axind] = int( np.ceil(oldshape[axind] / args['step']) )
        elif oplabel == 'thinBin':
            newshape[axind] = int( np.ceil(oldshape[axind] / args['bin']) )
        

        with h5py.File(dest.file, 'w') as df:
            destgrp = df[dest.group]
            
            #Create new data array
            destgrp.require_dataset('data', newshape, np.float32, chunks=True, compression='gzip')
            hdftools.copyAttrs(srcgrp['data'], destgrp['data'])
            
            #Copy the axes over, except the one being thinned
            for axis in dimlabels:
                if axis != ax:
                    srcgrp.copy(axis, destgrp)
                    
            #Create the new axis
            destgrp.require_dataset(ax, (newshape[axind],), np.float32, chunks=True)
            destgrp[ax].attrs['unit'] = srcgrp[ax].attrs['unit']
            
            
            #Initialize time-remaining printout
            #Chunks are big, so report more often than usual
            tr = util.timeRemaining(newshape[axind], reportevery=1)
   
            print(nchunks)
            for i in range(nchunks):
                print(i)
                #Update time remaining
                if verbose:
                    tr.updateTimeRemaining(i)
                sl= [ slice(None) ]*ndim
                
                if i != nchunks-1:
                    sl [chunkax] = slice(i*chunksize, (i+1)*chunksize, None)
                else:
                    sl [chunkax] = slice(i*chunksize, None, None)
                    
                
                if oplabel == 'avgDim':
                    avgDimOp(srcgrp['data'], destgrp['data'], sl, axind, args)
                elif oplabel == 'thinPick':
                    thinPickOp(srcgrp['data'], destgrp['data'], sl, axind, args)
                elif oplabel == 'thinBin':
                    pass
                else:
                    raise(ValueError, 'Unsupported oplabel! ' + str(oplabel))
                
                
            #Make the new axis
            if oplabel == 'avgDim':
                avgDimOp(srcgrp[ax], destgrp[ax], [slice(None)], 0, args)
            elif oplabel == 'thinPick':
                thinPickOp(srcgrp[ax], destgrp[ax], [slice(None)], 0, args)
            elif oplabel == 'thinBin':
                pass

            
            
            

    if delsrc:
        os.remove(src.file)
                









def thin(src, dest, ax, step=None, bin=None, loadwhole=False, 
         delsrc=False, verbose=False ):
    """
    'Thins' a dataset along a dimension with label 'ax'
    
    Parameters
    ----------
        src: hdfPath object
            Path to the source hdf dataset
            
        dest: hdfPath object
            Path to the destination hdf dataset

        ax: String
            Label of the axis to be thinned
            
        step: Int, default=10
            Thin the data by pulling out one datapoint every 'thin' indices.
            Thin is the default, and overrides bin.
            
        bin: Int
            Thin the dataset by averaging over bins of size 'bin'.
            
        loadwhole: Boolean
            If true, the entire hdf dataset will be loaded into mememory.
            If false, the data will be processed in chunks along the axis
            being thinned. Loading the whole file can be desirable if the file
            is relatively small and/or is chunked in an inconvenient way.
            
        delsrc: Boolean
            If this is set to true, the source file will be deleted after the
            run is complete.
            
        verbose: Boolean
            If this is set to true, printouts are enabled


    Returns
    -------
       dest: String
           Echos back the filepath to the newly created file if successful
    
    """
    
    #Coerece to type int
    if step is not None:
        step = int(step)
    if bin is not None:
        bin = int(bin)
    
    #Step overrides bin. Default is step=10
    if bin is None and step is None:
        step = 10
    elif bin is not None and step is not None:
        bin = None
    
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        
        validDataset(srcgrp)
        
        oldshape = list( srcgrp['data'].shape )
        ndim = len(oldshape)
        
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        
        try:
            axind = dimlabels.index(ax)
        except ValueError as e:
            print("ERROR: Couldn't find ax in dimensions array!: " + str(ax) )
            raise(e)

        #Calculate new shape
        newshape = np.copy(oldshape)
        
        if step is None:
            newshape[axind] = int( np.floor(oldshape[axind] / bin) )
        else:
            newshape[axind] = int( np.floor(oldshape[axind] / step) )


        with h5py.File(dest.file, 'w') as df:
            destgrp = df[dest.group]
            
            #Create new data array
            destgrp.require_dataset('data', newshape, np.float32, chunks=True, compression='gzip')
            hdftools.copyAttrs(srcgrp['data'], destgrp['data'])
            
            #Copy the axes over, except the one being thinned
            for axis in dimlabels:
                if axis != ax:
                    srcgrp.copy(axis, destgrp)
                    
            #Create the new axis
            destgrp.require_dataset(ax, (newshape[axind],), np.float32, chunks=True)
            destgrp[ax].attrs['unit'] = srcgrp[ax].attrs['unit']

             
            if step is not None and bin is None:
                #Thin by plucking out entries (faster, doesn't load whole array)
                b = int(np.floor(oldshape[axind]/step))*step
                srcslice = [ slice(None) ]*ndim
                srcslice[axind] = slice(0, b, step)
                srcaxslice = [srcslice[axind] ]
                
                if loadwhole:
                    srcdata = srcgrp['data']
                    destdata = srcdata[tuple(srcslice) ]
                    destgrp['data'][:] = destdata
                    destgrp[ax][:] = srcgrp[ax][ tuple(srcaxslice) ]  
                else:
                    destgrp['data'][:] = srcgrp['data'][ tuple(srcslice) ] 
                    destgrp[ax][:] = srcgrp[ax][ tuple(srcaxslice) ]
              
                
                
            elif step is None and bin is not None:
                
                if loadwhole:
                    srcdata = srcgrp['data']
                    destdata = np.empty(newshape)
                    destax = np.empty(newshape[axind])
                
                #Initialize time-remaining printout
                tr = util.timeRemaining(newshape[axind])
                
                for i in range(newshape[axind]):
                    
                    #Update time remaining
                    if verbose:
                        tr.updateTimeRemaining(i)
                    print(tr.avgExecutionTime())
                    
                    srcslice = [ slice(None) ]*ndim
                    srcslice[axind] = slice(i*bin, (i+1)*bin, None)
                    srcaxslice = [srcslice[axind] ]
                    
                    destslice = [ slice(None) ]*ndim
                    destslice[axind] = slice(i, i+1, None)
                    destaxslice = [ destslice[axind] ]
                    
                    #print(str(srcslice[axind]) + ' -> ' + str(destslice[axind]))
                    
                    
                    if loadwhole:
                        indrange = (i*bin, (i+1)*bin)
                        s = np.take(srcdata, indrange, axis=axind)
                        s = np.mean(s)
                        destdata[tuple(destslice)] = s
                    else:
                        s = srcgrp['data'][ tuple(srcslice) ] 
                        x = srcgrp[ax][tuple(srcaxslice)]
                        s = np.mean(s)
                        x = np.mean(x)
                        destgrp['data'][tuple(destslice)] = s
                        destgrp[ax][tuple(destaxslice)] = x
                        
                if loadwhole:
                    destgrp['data'][:] = destdata[:]
                    destgrp[ax][:] = destax[:]
            
    if delsrc:
        os.remove(src.file)
                

           
              
              

    
if __name__ == '__main__':
    full = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full.hdf5')
    thinned = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full_thinned.hdf5')
    avged = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full_avg.hdf5')
    
    x = thinPick(full, thinned, 'time', step = 10, verbose=True)
    #x = avgDim(full, avged, 'reps', verbose=True)