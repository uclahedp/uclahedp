#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 09:20:38 2019

@author: peter
"""
import h5py
import numpy as np
import os

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
    
    
    
def trimDim(src, dest, ax, bounds=None, values = None, delsrc=False, verbose=False):
    if bounds is None:
        print("Must set bounds in trimDim!")
        raise(ValueError)
    
    bounds = np.array(bounds, dtype=np.float32)
    
    if bounds[0] > bounds[1]:
        print("Flipping bounds so that the second element is larger!")
        bounds = np.reverse(bounds)
    
    if values is None:
        useind = True
    else:
        useind = False
        
    chunked_array_op(src, dest, ax, 'trimDim', delsrc=delsrc, verbose=verbose, 
                     bounds=bounds, useind=useind)
    
def trimDimOp(src_dset, dest_dset, sl, axind, args):
    #These have already been convered into indices at this point if needed
    a = args['bounds'][0]
    b = args['bounds'][1]

    sl[axind] = slice(a, b, None)
    s = src_dset[tuple(sl)]
    sl[axind] = slice(None, None, None)
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
    
    
    
def thinBin(src, dest, ax, bin=None, delsrc=False, verbose=False):
    if bin is None:
        bin = 10
    else:
        bin = int(bin)    
    chunked_array_op(src, dest, ax, 'thinBin', delsrc=delsrc, verbose=verbose, bin=bin)
    
    
def thinBinOp(src_dset, dest_dset, sl, axind, args):
    bin = args['bin']
    s = src_dset[tuple(sl)]
    nelm = s.shape[axind]
    nbins = int( np.ceil(nelm/bin) )
    
    dsl = np.copy(sl)
    for i in range(nbins):
        print('i=' + str(i))
        indrange = (i*bin, (i+1)*bin)
        x = np.take(s, indrange, axis=axind)
        x = np.mean(x, axis=axind)
        dsl[axind] = slice(None, None, None)
        dest_dset[tuple(dsl)] = x
        
    


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
        print("Chunking axis: " + str(dimlabels[chunkax]))
        
        if srcgrp['data'].chunks[chunkax] < 2:
            print("WARNING: POSSIBLE INEFFICENT CHUNKING DETECTED!")
        
            
        #Determine optimal chunksize (along chunkax)
        ideal_chunk_elms = 1e7 #1e7*4 bytes (per float32) ~ 40mb, which is good
        nper = np.product(oldshape)/oldshape[chunkax] #number of values per chunk ax value

        chunksize = int( np.round(ideal_chunk_elms/nper) )
        if chunksize < 1:
            chunksize = 1
        
        #Determine nchunks    
        nchunks = int( np.ceil(oldshape[chunkax] / chunksize))
        

        #Make alterations to the shape for the new dataset
        #This will obviously depend on what operation is being performed
        newshape = np.copy(oldshape)
        if oplabel == 'avgDim':
            newshape[axind] = 1
            opfunc = avgDimOp
        elif oplabel == 'thinPick':
            newshape[axind] = int( np.ceil(oldshape[axind] / args['step']) )
            opfunc = thinPickOp
        elif oplabel == 'thinBin':
            newshape[axind] = int( np.ceil(oldshape[axind] / args['bin']) )
            opfunc = thinBinOp
        elif oplabel == 'trimDim':
            #If values are being passed, figure out the indices here
            if not args['useind']:
                if args['bounds'][0] < srcgrp[ax][:].min():
                    a = 0
                else:
                    a  =  np.abs(srcgrp[ax][:] - args['bounds'][0]).argmin()
                
                if args['bounds'][1] > srcgrp[ax][:].max():
                    b = oldshape[axind] -1
                else:
                    b  =  np.abs(srcgrp[ax][:] - args['bounds'][1]).argmin() 
                args['bounds'] = (a, b)
            
            args['bounds'] = np.clip(args['bounds'], 0, oldshape[axind]-1)
            newshape[axind] = np.abs(args['bounds'][1] - args['bounds'][0] )
            opfunc = trimDimOp
        else:
            raise(ValueError, 'Unsupported oplabel! ' + str(oplabel))
        
        with h5py.File(dest.file, 'w') as df:
            destgrp = df[dest.group]
			
			#Copy all the dataset attributes
			hdftools.copyAttrs(srcgrp, destgrp)
            
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
            tr = util.timeRemaining(nchunks, reportevery=1)

            for i in range(nchunks):
                #Update time remaining
                if verbose:
                    tr.updateTimeRemaining(i)
                sl= [ slice(None) ]*ndim
                
                if i != nchunks-1:
                    sl [chunkax] = slice(i*chunksize, (i+1)*chunksize, None)
                else:
                    sl [chunkax] = slice(i*chunksize, None, None)
                    
               
                opfunc(srcgrp['data'], destgrp['data'], sl, axind, args)

            #Make the new axis
            opfunc(srcgrp[ax], destgrp[ax], [slice(None)], 0, args)
     
    if delsrc:
        os.remove(src.file)
                


    
if __name__ == '__main__':
    
    #Windows
    full = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full.hdf5") )
    thinned = hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full_thinned.hdf5") )
    trimmed = hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full_trim.hdf5") )
    avged = hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full_avg.hdf5") )
    
    #OSX
    #full = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full.hdf5')
    #thinned = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full_thinned.hdf5')
    #avged = hdftools.hdfPath('/Volumes/PVH_DATA/LAPD_Mar2018/FULL/run61_LAPD1_full_avg.hdf5')
    
    x = trimDim(full, trimmed, 'time', bounds=[0,1e-5], values=False, verbose=True)
    #x = thinBin(full, thinned, 'time', bin = 10, verbose=True)
    #x = thinPick(full, thinned, 'time', step = 10, verbose=True)
    #x = avgDim(full, avged, 'reps', verbose=True)