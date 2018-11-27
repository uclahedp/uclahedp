#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rawtools.py: Tools for creating HDF5 files formatted in the UCLA HEDP "Raw" manner

Created by Scott Feister on Tue Nov 20 08:22:44 2018
"""

from scipy.io import readsav
import os
import h5py
import numpy as np
import datetime

def sraw2hraw(fname_sav):
    """ Convert an IDL-Sav Raw file into an HDF5 Raw file

    Parameters
    ----------
        fname_sav : str
            Name of the IDL-Sav Raw file

    Returns
    -------
        fname_h5
            Name of the HDF5 Raw file
    """
    
    # Read the IDL raw save file into a python dictionary, d
    idl_dict = readsav(fname_sav, python_dict=True)
    mynames = idl_dict['struct'].dtype.names
    myitems = idl_dict['struct'][0]
    
    d = {} # Temporary dictionary holding values
    for name, item in zip(mynames, myitems): # Python3 syntax
        if isinstance(item, bytes): # Interpret any byte strings as UTF-8
            item = item.decode('UTF-8')
        d[name] = item
    
    # Homogenize instances of "None" for case consistency
    for k in d.keys():
        if isinstance(d[k], str) and (d[k] == 'None' or d[k] == 'none'):
            d[k] = str(None)
    
    # Create the HDF5 raw save filename
    fname_h5 = os.path.splitext(fname_sav)[0] + ".h5"
    
    # Open the HDF5 raw save file and write IDL elements per the new syntax
    with h5py.File(fname_h5, 'w') as f:
        # Specify the HDF5 file creation time
        f.attrs['H5_creation_time'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.attrs['IDL_raw_filename'] = str(fname_sav)
        
        # Specify the probe parameters
        f.attrs['gridded'] = True
        f.attrs['probe_name'] = d['PROBE']
        f.attrs['probe_type'] = str(None)

        goodkeys = ('DAQ', 'DRIVE', 'PLANE', 'DATA_TYPE', 'DATA_FORM', 'RUN', 'NOTES')
        for k in goodkeys: # Add everything else into the HDF5 file as an attribute
            f.attrs[k.lower()] = d[k]

        f.attrs['dt_secs'] = d['DT']
        f.attrs['time_unit'] = 'us' # TODO
        f.attrs['space_unit'] = 'cm' # TODO
        f.attrs['chan_labels'] = str(None) # TODO
        f.attrs['pos_labels'] = [a.encode('utf-8') for a in ['X', 'Y', 'Z']] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        
        # Extract info needed to create "data" and "pos" arrays
        idl_data = d['DATA'].T # Transpose puts it back into correct index order
        
        [nti, nx, ny, nz, nreps, nchan] = idl_data.shape
        npos = nx * ny * nz # Total number of positions recorded
        f['data'] = np.reshape(idl_data, [nti, npos, nreps, nchan])
        ndim = 3
        pos = np.zeros([ndim, npos])
        xgv = d['XAXES']
        ygv = d['YAXES']
        zgv = d['ZAXES']
        X, Y, Z = np.meshgrid(xgv, ygv, zgv, indexing='ij')
        pos[0,:] = X.flatten()
        pos[1,:] = Y.flatten()
        pos[2,:] = Z.flatten()
        f['pos'] = pos

    return d

if __name__ == "__main__":
    # Quick and dirty test of the array manipulations
    idl_data = np.array([[0,2,3],[5,1,2]])
    xgv = np.array(['19', '20'])
    ygv = np.array(['30', '40', '50'])
    (nx, ny) = idl_data.shape
    data =  np.reshape(idl_data, [nx*ny])
    X, Y = np.meshgrid(xgv, ygv, indexing='ij')
    pos = np.zeros([2, nx*ny])
    pos[0,:] = X.flatten()
    pos[1,:] = Y.flatten()

    fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_LAPD1_pos_raw.sav"
    #fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_tdiode_t_raw.sav"
    newdict = sraw2hraw(fname_sav)

