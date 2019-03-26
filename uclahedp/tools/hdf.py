#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import h5py
import os
import numpy as np

from uclahedp import csvtools, hdftools







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
        attrs[k] = ( csvtools.fixType(obj.attrs[k][0]), obj.attrs[k][1])
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
    
    src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_LAPD1.hdf5")  )
    with h5py.File(src.file, 'r') as f:
        datagrp = f[src.group]
        print( validDataset(datagrp) )
