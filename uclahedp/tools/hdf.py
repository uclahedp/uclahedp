#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import h5py
import os
import numpy as np

import pathlib

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools



def copyDataset(sf, df, name):
    num = sf[name].shape
    dtype = sf[name].dtype
    df.require_dataset(name, num, dtype, chunks=True)[:] = sf[name][:]
    copyAttrs(sf[name], df[name])

def requireDirs(destfile):
    """
    Checks to see if the folder associated with a filepath exists, and creates
    it if it does not.
    """
    p = pathlib.Path(destfile)
    folder = p.parent
    if not folder.exists():
        os.mkdir(folder)


def writeAttrs(attrs, group):
    """
    Writes a dictionary of attributes (attrs) onto the attrs of an HDF5
    Group (group). Any conflicts will be overwritten.
    
    This function assumes the attributes to be tuples of (value, unit) where
    value may be an int,float, or str, and unit is a str. 
    """
    for k in attrs:
            if attrs[k] is not None:
                val, unit =  attrs[k]
                val = str(val).encode('utf-8')
                unit = str(unit).encode('utf-8')
                group.attrs.create(k,  (val, unit) )
                
def readAttrs(obj):
    """
    Reads attributes from an HDF5 object and returns them as a dictionary
    """
    attrs = {}
    for k in obj.attrs.keys():
        val = obj.attrs[k][0].decode('utf-8')
        unit = obj.attrs[k][1].decode('utf-8')
        

        val = csvtools.fixType(val)
        
        attrs[k] = (val,unit)
    return attrs


def copyAttrs(srcobj, destobj):
    """
    Copies attributes from one HDF5 object to another, overwriting any
    conflicts
    """
    for k in srcobj.attrs.keys():
                destobj.attrs[k] = srcobj.attrs[k]
    return True


def arrToStrList(arr):
    l = arr.tolist()
    return [x.decode('utf-8') for x in l]
    
def strListToArr(strlist):
    return np.array( [s.encode('utf-8') for s in strlist] )


class hdfPath():
    """
    The hdfPath class simply allows convenient passing of a path to both an
    HDF file and a group within that file.
    
    The __str__ method enables the path to be easily printed as a str.
    """
    def __init__(self, file, group=''):
        self.file = file
        self.group = group + '/'
        
    def __str__(self):
        return str(self.file) + ' : ' + str(self.group)
    
    
    




class hdfDatasetExists(Exception):
    """A dataset by this name already exists in this group.\n 
        Due to design limitations of the HDF5 format, deleting/overwriting datasets is not recommended. \n
        Please either rename the dataset, or delete the entire hdf file and try again.\n
    """
    pass

class hdfGroupExists(Exception):
    """A group by this name already exists.\n 
        Due to design limitations of the HDF5 format, deleting/overwriting groups is not recommended. \n
        Please either rename the group, or delete the entire hdf file and try again.\n
    """
    pass





if __name__ == "__main__":
    
    exp = 'LAPD_Jan2019'
    probe = 'LAPD10'
    run = 29
    
    src = hdftools.hdfPath( '/Volumes/PVH_DATA/' + exp + '/RAW/run' + str(run) + '_' + probe + '_raw.hdf5')
    
    print(src.file)
    with h5py.File(src.file, 'r') as f:
        srcgrp = f[src.group]
        attrs = readAttrs(srcgrp)
        print(len(attrs))
