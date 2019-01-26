#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import h5py
import csvtools

def valid_dataset(dataset):
    
    ndims = len( [d for d in dataset.dims ])
    shape = dataset.shape

    if len(shape ) != ndims :
        raise ValueError("Wrong number of dimensions!")

        
    for i,d in enumerate(dataset.dims):
        #print(type(d[0]))
        label = d.label
        
        print(d[0])
        
        #print(shape[i])
        #print(d[0].shape[0])
        if shape[i] != d[0].shape[0]:
            raise ValueError("Dimensionsionality conflict")
    
    
    return True



def writeAttrs(attrs, group):
    for k in attrs:
            if attrs[k] is not None:
                val, unit =  attrs[k]
                val = str(val).encode('utf-8')
                unit = str(unit).encode('utf-8')
                group.attrs.create(k,  (val, unit) )
                
def readAttrs(obj):
    attrs = {}
    for k in obj.attrs.keys():
        attrs[k] = ( csvtools.fixType(obj.attrs[k][0]), obj.attrs[k][1])
    return attrs


def copyAttrs(srcobj, destobj):
    for k in srcobj.attrs.keys():
                destobj.attrs[k] = srcobj.attrs[k]
    return True




class hdfPath():
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
    
    filepath = '/Users/peter/Desktop/testhdf.h5'
    dataset_path = '/data'
    with h5py.File(filepath, 'r') as f:
        dataset = f['/data']
        print(valid_dataset(dataset))