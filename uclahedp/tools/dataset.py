#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 09:20:38 2019
@author: peter

The dataset tools package contains functions for manupulating datasets that
follow the UCLA HEDP dataset format (defined in the validDataset function).

**Chunked Operations**
Many operations are wrapped in the chunked_array_op function to efficiently
chunk them to minimize memory usage. Each of these operations uses three
function levels. For example, avgDim:
    
avgDimOp -> This function does the actual operation (averaging) on a slice
of a dataset.

chunked_array_op -> This function determines the optimal chunking and
repeatedly calls avgDimOp on slices of the data until the entire operation is
completed. 

avgDim -> The user-level function. Takes some parameters, then calls
chunked_array_op to handle the operation

**delscr**
Overwriting data inside an HDF5 file is not supported because of limitations
with the HDF format: deleted data does not free memory, so file size ballons.
As an alternative, these operations have a delsrc keyword that will delete the
source dataset once the operation is complete. The operations can then be
'chained' together to avoid accumulating many copies of the data on disk.

"""
import h5py
import numpy as np
import os

from uclahedp.tools import hdf as hdftools
from uclahedp.tools import util

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


def getAxInd(ax, dimlabels):
    #Get ax index
    try:
        axind = dimlabels.index(ax)
    except ValueError as e:
        print("ERROR: Couldn't find ax in dimensions array!: " + str(ax) )
        raise(e)
    return axind



def avgDim(src, dest, ax, delsrc=False, verbose=False):
    """
    Average over one dimension of a dataset (collapsing it to len=1)
    src -> Source dataset (hdfpath object)
    dest -> Destination dataset path (hdfpath object)
    ax -> Axis to apply op to (name)
    delsrc -> Boolean, if true src file will be deleted after operation
    verbose -> Boolean, if true activates printouts
    """
    
    #Load some file parameters to calculate the shape of the new dataset   
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        oldshape = srcgrp['data'].shape
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        #Get ax index
        axind = getAxInd(ax, dimlabels)
        newshape = np.copy(oldshape)
        newshape[axind] = 1
    
    #Call the avgDim function, wrapped in the chunked_array_op framework
    chunked_array_op(src, dest, ax, avgDimOp, newshape, delsrc=delsrc, verbose=verbose)
    
def avgDimOp(src_dset, dest_dset, sl, axind, args):
    """
    Average over one dimension of a dataset (collapsing it to len=1)
    src_dset -> Source dataset (hdfpath to 'data')
    dest_dset -> Destination dataset path (hdfpath to 'data')
    sl -> Slice of source data to apply current operation to (for chunks)
    axind -> Axis to apply operation to
    args -> Additional function args passed from higher function
    """
    s = src_dset[tuple(sl)]
    s = np.mean(s, axis=axind, keepdims=True)
    sl[axind] = slice(None, None, None)  # This dim is 1 in the dest dataset
    dest_dset[tuple(sl)] = s


def trimDim(src, dest, ax, ind_bounds=None, val_bounds = None, delsrc=False, verbose=False):
    """
    Trim a dimension of a dataset, disgarding some data
    
    src -> Source dataset (hdfpath object)
    dest -> Destination dataset path (hdfpath object)
    ax -> Axis to apply op to (name)
    bounds -> Start and stop bounds for trim. Default is indicies, but
    interpreted as values if 'values' flag is set.
    values -> Boolean, if true interpret bounds as axis values not indices.
    delsrc -> Boolean, if true src file will be deleted after operation
    verbose -> Boolean, if true activates printouts
    """

        
    if not ind_bounds and not val_bounds:
        print("Using ind bounds (as default)")
        bounds = (None, None)
    if ind_bounds and not val_bounds:
        print("Using ind bounds")
        bounds = ind_bounds
    elif val_bounds and not ind_bounds:
        print("Using val bounds")
        #If values are being passed, figure out the indices here
        with h5py.File(src.file, 'r') as sf:
            srcgrp = sf[src.group]
            oldshape = srcgrp['data'].shape
            dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
            #Get ax index
            axind = getAxInd(ax, dimlabels)
            
            if val_bounds[0] < srcgrp[ax][:].min():
                a = 0
            else:
                a  =  np.abs(srcgrp[ax][:] - val_bounds[0]).argmin()
                
            if val_bounds[1] > srcgrp[ax][:].max():
                b = oldshape[axind] -1
            else:
                b  =  np.abs(srcgrp[ax][:] - val_bounds[1]).argmin() 
            bounds = (a, b)
            bounds = np.clip(bounds, 0, oldshape[axind]-1)
            print(bounds)
    else:
        raise ValueError("Cannot specify ind_bounds AND val_bounds!")
        
        
    #Load some file parameters to calculate the shape of the new dataset   
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        oldshape = srcgrp['data'].shape
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        #Get ax index
        axind = getAxInd(ax, dimlabels)
        newshape = np.copy(oldshape)
        newshape[axind] = np.abs(bounds[1] - bounds[0] )
        
        
        
        
    chunked_array_op(src, dest, ax, trimDimOp, newshape, delsrc=delsrc, verbose=verbose, 
                     bounds=bounds)
    
def trimDimOp(src_dset, dest_dset, sl, axind, args):
    """
    Trim a dimension of a dataset, disgarding some data
    
    src_dset -> Source dataset (hdfpath to 'data')
    dest_dset -> Destination dataset path (hdfpath to 'data')
    sl -> Slice of source data to apply current operation to (for chunks)
    axind -> Axis to apply operation to
    args -> Additional function args passed from higher function
    """
    
    #These have already been convered into indices at this point if needed
    a = args['bounds'][0]
    b = args['bounds'][1]

    sl[axind] = slice(a, b, None)
    s = src_dset[tuple(sl)]
    sl[axind] = slice(None, None, None)
    dest_dset[tuple(sl)] = s
    

#Thin by picking out values (not by binning and averaging)
def thinPick(src, dest, ax, step=None, delsrc=False, verbose=False):
    """
    Thin a dataset by picking every nth point and disgarding the rest
    
    src -> Source dataset (hdfpath object)
    dest -> Destination dataset path (hdfpath object)
    ax -> Axis to apply op to (name)
    step -> The points kept will be indices i*step
    delsrc -> Boolean, if true src file will be deleted after operation
    verbose -> Boolean, if true activates printouts
    """
    if step is None:
        step = 10
    else:
        step = int(step)
        
    #Load some file parameters to calculate the shape of the new dataset   
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        oldshape = srcgrp['data'].shape
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        #Get ax index
        axind = getAxInd(ax, dimlabels)
        newshape = np.copy(oldshape)
        newshape[axind] = int( np.ceil(oldshape[axind] / step) )
        
    chunked_array_op(src, dest, ax, thinPickOp, newshape, delsrc=delsrc, verbose=verbose, step=step)

def thinPickOp(src_dset, dest_dset, sl, axind, args):
    """
    Thin a dataset by picking every nth point and disgarding the rest
    
    src_dset -> Source dataset (hdfpath to 'data')
    dest_dset -> Destination dataset path (hdfpath to 'data')
    sl -> Slice of source data to apply current operation to (for chunks)
    axind -> Axis to apply operation to
    args -> Additional function args passed from higher function
    """
    #Thin by plucking out entries
    sl[axind] = slice(None, None, args['step'])
    s = src_dset[tuple(sl)]
    sl[axind] = slice(None, None, None)
    print(s.shape)
    dest_dset[tuple(sl)] = s
    
    
#Thin by binning and averaging
def thinBin(src, dest, ax, bin=None, delsrc=False, verbose=False):
    """
    Thin a dataset by averaging it over non-overlapping bins.
    
    src -> Source dataset (hdfpath object)
    dest -> Destination dataset path (hdfpath object)
    ax -> Axis to apply op to (name)
    bin -> The width of each bin
    delsrc -> Boolean, if true src file will be deleted after operation
    verbose -> Boolean, if true activates printouts
    """
    
    if bin is None:
        bin = 10
    else:
        bin = int(bin)
        
    #Load some file parameters to calculate the shape of the new dataset   
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        oldshape = srcgrp['data'].shape
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        #Get ax index
        axind = getAxInd(ax, dimlabels)
        newshape = np.copy(oldshape)
        newshape[axind] = int( np.ceil(oldshape[axind] / bin) )
    
    chunked_array_op(src, dest, ax, thinBinOp, newshape, delsrc=delsrc, verbose=verbose, bin=bin)
    
    
def thinBinOp(src_dset, dest_dset, sl, axind, args):
    """
    Thin a dataset by averaging it over non-overlapping bins.
    
    src_dset -> Source dataset (hdfpath to 'data')
    dest_dset -> Destination dataset path (hdfpath to 'data')
    sl -> Slice of source data to apply current operation to (for chunks)
    axind -> Axis to apply operation to
    args -> Additional function args passed from higher function
    """
    
    bin = args['bin']
    s = src_dset[tuple(sl)]
    nelm = s.shape[axind]
    nbins = int( np.ceil(nelm/bin) )
    
    dsl = np.copy(sl)

    
    for i in range(nbins):
        indrange = np.arange(i*bin, (i+1)*bin)
        x = np.take(s, indrange, axis=axind)
        x = np.mean(x, axis=axind)
        #Should this maybe be (i, i+1, None)?
        dsl[axind] = slice(i, i+1, None)
        dest_dset[tuple(dsl)] = x
        
    


def chunked_array_op(src, dest, ax, op, newshape, delsrc=False, verbose=False, **args):
    """
    Apply one of the array functions to an entire dataset, breaking the
    dataset up into chunks to keep memory load low.
    
    src -> Source dataset (hdfpath object)
    dest -> Destination dataset path (hdfpath object)
    ax -> Axis (0 indexed) to average
    op -> Function to be applied. This function must be one of the op functions
    defined in this file, and must be included in the elif tree in this function
    newshape -> Shape the new dataset will be after op has been applied
    delsrc -> Boolean, if true src file will be deleted after operation
    verbose -> Boolean, if true activates printouts
    """
    
    with h5py.File(src.file, 'r') as sf:
        srcgrp = sf[src.group]
        
        #Check source is valid dataset
        validDataset(srcgrp)
        
        #Load information about source dataset
        oldshape = list( srcgrp['data'].shape )
        ndim = len(oldshape)
        dimlabels = hdftools.arrToStrList( srcgrp['data'].attrs['dimensions'][:] )
        
        #Get ax index
        axind = getAxInd(ax, dimlabels)
        
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
        

        #Create the destination dataset
        with h5py.File(dest.file, 'w') as df:
            destgrp = df[dest.group]
            
            #Copy all the dataset attributes
            hdftools.copyAttrs(srcgrp, destgrp)
            
            #Create new data array
            destgrp.require_dataset('data', newshape, np.float32, chunks=True, compression='gzip')
            hdftools.copyAttrs(srcgrp['data'], destgrp['data'])
            
            #Copy the axes over, except the one being operated on
            #That axis will be copied over later, with changes
            for axis in dimlabels:
                if axis != ax:
                    srcgrp.copy(axis, destgrp)
                    
            #Create the axis being operated on
            #Newshape was determined above, and is specific to the op
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
                
                #Assemble the chunk slices
                if i != nchunks-1:
                    sl [chunkax] = slice(i*chunksize, (i+1)*chunksize, None)
                else:
                    sl [chunkax] = slice(i*chunksize, None, None)
                    
               #Apply op to the chunk
                op(srcgrp['data'], destgrp['data'], sl, axind, args)

            #Make the new axis by applying op to the old axis
            op(srcgrp[ax], destgrp[ax], [slice(None)], 0, args)
     
    #If requested, delete the source file
    if delsrc:
        os.remove(src.file)
                


    
if __name__ == '__main__':
    
    #Windows
    #full = hdftools.hdfPath( os.path.join("F:", os.sep, "2019BIERMANN","FULL", "run29_LAPD_C6_full.hdf5") )
    #thinned = hdftools.hdfPath(os.path.join("F:", os.sep, "2019BIERMANN","FULL", "run29_LAPD_C6_avg.hdf5") )
    #trimmed = hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full_trim.hdf5") )
    #avged = hdftools.hdfPath(os.path.join("F:", "LAPD_Mar2018", "FULL", "run61_LAPD1_full_avg.hdf5") )
    
    #OSX
    full = hdftools.hdfPath('/Volumes/PVH_DATA/2019BIERMANN/FULL/run29_LAPD_C6_full.hdf5')
    thinned = hdftools.hdfPath('/Volumes/PVH_DATA/2019BIERMANN/FULL/run29_LAPD_C6_avg.hdf5')
    trimmed = hdftools.hdfPath('/Volumes/PVH_DATA/2019BIERMANN/FULL/run29_LAPD_C6_trim.hdf5')
    avged = hdftools.hdfPath('/Volumes/PVH_DATA/2019BIERMANN/FULL/run29_LAPD_C6_dimavged.hdf5')
    
    #x = trimDim(full, trimmed, 'time', val_bounds=[0,1e-6],verbose=True)
    #x = trimDim(full, trimmed, 'time', ind_bounds=[120, 500],verbose=True)
    
    x = thinBin(full, thinned, 'shots', bin = 5, verbose=True)
    #x = thinPick(full, thinned, 'shots', step = 10, verbose=True)
    #x = avgDim(full, avged, 'shots', verbose=True)