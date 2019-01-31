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

# TODO: This function should grid the t0 vector too? Maybe use same function but with a keyword eg. quantity='data'?
def regrid(f):
    """ Put the data in the open HDF5 file "f" back onto a grid
    Parameters
    ----------
        f : file object
            Opened file object (allowing read mode) for an HDF5 Raw file       
    Returns
    -------
        gridded : numpy.ndarray
            Gridded data of shape [nti, nx, ny, nz, npos, nreps, nchan]
        xgv, ygv, zgv : numpy.ndarray
            1D axis vectors for x, y, and z
    """

    if not f.attrs['gridded']:
        raise Exception("Data not gridded!")

    nx = f.attrs['grid_nx']
    ny = f.attrs['grid_ny']
    nz = f.attrs['grid_nz']

    [nti, npos, nreps, nchan] = f['data'].shape
    gridded = np.reshape(f['data'], [nti, nx, ny, nz, nreps, nchan])

    xgv = f['grid_xgv']
    ygv = f['grid_ygv']
    zgv = f['grid_zgv']

    return gridded, xgv, ygv, zgv


#TODO: move this time vector function into a file of functions that apply to ANY HEDP HDF5 file, not just raw ones?
#TODO: handle lists of positions, reps and produce an array of time vectors all in one function call?
def timevector(f, pos=0, rep=0):
    """ Generates a time vector for a given shot in any HDF5 file
    Parameters
    ----------
        f : file object
            Opened file object (allowing read mode) for an HDF5 Raw file

        pos : integer
            Position index of the requested time vector

        rep : integer
            Repetition index of the requested time vector
    Returns
    -------
       time : Time vector in time units of the base HDF file, shifted by
       the appropriate time shift.
    """
    time = np.arange(f.attrs['nti'])*f.attrs['dt'] - f['t0'][pos, rep]

    return time


def sraw2hraw(fname_sav):
    """ Convert an IDL-Sav Raw file into an HDF5 Raw file
    Accepts data_forms of both "t" (diode reading) and "pos" (position scan)
    For "pos" types, requires the items "XAXES", "YAXES", and "ZAXES"
    in the IDL struct

    Parameters
    ----------
        fname_sav : str
            Name of the IDL-Sav Raw file

    Returns
    -------
        fname_h5 : str
            Name of the HDF5 Raw file
    """

    # Read the IDL raw save file into a python dictionary, d
    idl_dict = readsav(fname_sav, python_dict=True)
    mynames = idl_dict['struct'].dtype.names
    myitems = idl_dict['struct'][0]

    d = {}  # Temporary dictionary holding values
    for name, item in zip(mynames, myitems):  # Python3 syntax
        if isinstance(item, bytes):  # Interpret any byte strings as UTF-8
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
        
        # CHANGE: Removed these two attributes, replace them with a log list that can be appended
        # each time a file is changed, maitaining a paper trail of what functions have touched it.
        # Specify HDF5 file creation parameters
        #f.attrs['H5_creation_time'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        #f.attrs['IDL_raw_filename'] = str(fname_sav)
        
        entry = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': Created from IDL sav file: ' + str(fname_sav)
        log = [entry]
        f['log'] = [s.encode('utf-8') for s in log] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289

        # Specify probe parameters
        f.attrs['probe_name'] = d['PROBE']
        f.attrs['probe_type'] = str(None)  # Not included in the IDL save files
        
        # Specify some time parameters
        f.attrs['dt'] = d['DT']  #Store in the same units as the t0 time for consistency
        f.attrs['dt_unit'] = 's'
        f.attrs['chan_labels'] = [s for s in d['CHAN_TITLE']]  # Note that these elements are already in utf-8 format

        # Specify some spatial parameters
        if d['DATA_FORM'] == "t":
            f.attrs['space_unit'] = str(None)
            f.attrs['pos_labels'] = str(None)
            f.attrs['gridded'] = False
        else:
            f.attrs['space_unit'] = 'cm'  # TODO
            f.attrs['pos_labels'] = [s.encode('utf-8') for s in ['X', 'Y', 'Z']] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
            f.attrs['gridded'] = True

        # Specify several additional parameters
        goodkeys = ('DAQ', 'DRIVE', 'PLANE', 'DATA_TYPE',
                    'DATA_FORM', 'RUN', 'NOTES')
        for k in goodkeys: # Add everything else into the HDF5 file as an attribute
            if k in d.keys():
                f.attrs[k.lower()] = d[k]
            else:
                f.attrs[k.lower()] = str(None)

        # Extract info needed to create "data" and "pos" arrays
        idl_data = d['DATA'].T # Transpose puts it back into correct index order

        # Write the "pos" array if appicable, obtain array dims
        if d['DATA_FORM'] == "t":
            [nti, nreps, nchan] = idl_data.shape
            npos = 1
        else:
            ndim = 3
            [nti, nx, ny, nz, nreps, nchan] = idl_data.shape
            npos = nx * ny * nz  # Total number of positions recorded
            f.attrs['npos'] = npos
            pos = np.zeros([ndim, npos])
            xgv = d['XAXES']
            ygv = d['YAXES']
            zgv = d['ZAXES']
            X, Y, Z = np.meshgrid(xgv, ygv, zgv, indexing='ij')
            pos[0, :] = X.flatten()
            pos[1, :] = Y.flatten()
            pos[2, :] = Z.flatten()
            #f['grid_xgv'] = xgv
            #f['grid_xgv'].attrs['unit'] = 'cm'
            #f['grid_ygv'] = ygv
            #f['grid_ygv'].attrs['unit'] = 'cm'
            #f['grid_zgv'] = zgv
            #f['grid_zgv'].attrs['unit'] = 'cm'
            #f['pos'] = pos
            #f['pos'].attrs['unit'] = 'cm'
            #f.attrs['grid_nx'] = nx
            #f.attrs['grid_ny'] = ny
            #f.attrs['grid_nz'] = nz
            
            
           
            
            
        

        # Save the number of timesteps, reps, and channels for later
        # easy access.
        #f.attrs['nti'] = nti
        #f.attrs['nreps'] = nreps
        #f.attrs['nchan'] = nchan

        # Write the "data" array
        f['data'] =  np.squeeze( np.reshape(idl_data, [nti, nx, ny, nz, nreps, nchan]) )
        #f['data'].attrs['unit'] = 'V'
        # Create the t0 array
        # This array stores the time
        #f['t0'] = np.zeros([npos, nreps])
        #f['t0'].attrs['unit'] = 's'
        

        # These parts of the hdf file are dependent on what the dimensions are
        dimlabels = []
        units=[]
        axnames = ['ax0', 'ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6', 'ax7', 'ax8']
        if nti is not 1:
            dimlabels.append('t')
            units.append('s')
            f[axnames.pop(0)] = np.arange(nti)*f.attrs['dt']
            
        if nx is not 1:
            dimlabels.append('x')
            units.append('cm')
            f[axnames.pop(0)] = xgv
        else:
            f.attrs['x0'] = xgv[0]
            
        if ny is not 1:
            dimlabels.append('y')
            units.append('cm')
            f[axnames.pop(0)] = ygv
        else:
            f.attrs['y0'] = ygv[0]
            
        if nz is not 1:
            dimlabels.append('z')
            units.append('cm')
            f[axnames.pop(0)] = zgv
        else:
            f.attrs['z0'] = zgv[0]

        if nreps is not 1:
            dimlabels.append('reps') 
            units.append('') # Astropy.unit will interpret empty string as dimensionless
            f[axnames.pop(0)] = np.arange(nreps)
        if nchan is not 1:
            dimlabels.append('channels')
            units.append('')
            f[axnames.pop(0)] = np.arange(nchan)
            
        f.attrs['dimlabels'] = [s.encode('utf-8') for s in dimlabels] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        f.attrs['dimunits'] = [s.encode('utf-8') for s in units] # Note 'utf-8' syntax is a workaround for h5py issue: https://github.com/h5py/h5py/issues/289
        f.attrs['data_label'] = 'Raw Data'
        f.attrs['data_unit'] = 'V'

        return fname_h5


if __name__ == "__main__":

    #fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_LAPD1_pos_raw.sav"
    #fname_sav = r"C:\Users\scott\Documents\DATA\2018-11-26 Example UCLA Raw files\run56_LAPD1_pos_raw.sav"
    #fname_sav = r"C:\Users\scott\Documents\UCLA\IDL to Python Bdot\DataForScott\DataForScott\RAW\run40_tdiode_t_raw.sav"
    #fname_sav = r"C:\Users\scott\Documents\DATA\2018-11-26 Example UCLA Raw files\run102_PL11B_pos_raw.sav"
    fname_sav = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.sav"
    
    #fname_sav = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run102_PL11B_pos_raw.sav"
    
    
    fname_sav = r"C:\Users\Peter\Desktop\TempData\run102_PL11B_pos_raw.sav"
    
    fname_h5 = sraw2hraw(fname_sav)

    #with h5py.File(fname_h5, 'r') as f:
    #    gridded, xgv, ygv, zgv = regrid(f)
    #    print(np.shape(gridded))
        
    print('Done')

    
    

