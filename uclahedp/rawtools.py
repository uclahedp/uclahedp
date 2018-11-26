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
    
    idl_dict = readsav(fname_sav, python_dict=True)
    mynames = idl_dict['struct'].dtype.names
    myitems = idl_dict['struct'][0]
    
    d = {} # Temporary dictionary holding values
    for name, item in zip(mynames, myitems): # Python3 syntax
        if isinstance(item, bytes): # Interpret any byte strings as UTF-8
            item = item.decode('UTF-8')
        d[name] = item
    
    
    
    fname_h5 = os.path.splitext(fname_sav)[0] + ".h5"
    
    print(d)
    with h5py.File(fname_h5, 'w') as f:
        # Specify the HDF5 file creation time
        f.attrs["HDF5CreationTime"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Extract info needed to create "data" array
        nti = d['NTI']
        nx = d['NXPOS']
        ny = d['NYPOS']
        nz = d['NZPOS']
        nreps = d['NREPS']
        nchan = d['NCHAN']
        
        # Allocate the "data" array as zeros
        f['data'] = np.zeros((nti, nx, ny, nz, nreps, nchan), 'f')
        
        # Specify axes grid vectors (except reps and chans)
        f["tgv"] = d['TIME']
        f["xgv"] = d['XAXES']
        f["ygv"] = d['YAXES']
        f["zgv"] = d['ZAXES']
        
        # Specify axes labels
#        f.attrs["tlabel"] = "Time"
#        f.attrs["xlabel"] = "X position"
#        f.attrs["ylabel"] = "Y position"
#        f.attrs["zlabel"] = "Z position"       

        f.attrs["tlabel"] = d['TIME_TITLE']
        f.attrs["xlabel"] = d['AXES_TITLE'][0].decode('UTF-8')
        f.attrs["ylabel"] = d['AXES_TITLE'][1].decode('UTF-8')
        f.attrs["zlabel"] = d['AXES_TITLE'][2].decode('UTF-8')
        f.attrs["repslabel"] = "Repetition number"
        f.attrs["chanlabel"] = "Channel"

        usedkeys = ('DATA', 'TIME', 'XAXES', 'YAXES', 'ZAXES', 'CHAN_TITLE', 'TIME_TITLE', 'AXES_TITLE')
        for k in (d.keys() - usedkeys): # Add everything else into the HDF5 file as an attribute
            f.attrs[k] = d[k]

        
    #    for k, v in newdict.items():
    #        print(k)
    #        f[k] = v
    #for item in myitems:
    #    print(type(item))
    
    return d

if __name__ == "__main__":
    fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_LAPD1_pos_raw.sav"
    #fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_tdiode_t_raw.sav"
    newdict = sraw2hraw(fname_sav)
