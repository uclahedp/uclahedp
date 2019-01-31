#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: peter
Example 01
Demonstrates use of the lapdtools and bdot programs to load in and process
an LAPD dataset. The metadata CSVs are included.

THE HDF FILE FOR THIS EXAMPLE IS NOT INCLUDED DUE TO SIZE

This dataset is a single Raptor shot from the Mar2018 LAPD experiment.
"""

import os
import hdftools
import lapdtools
import bdot

#Get the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))

#Define paths to two directories
#If you leave the script in the path it downloads in, these should already
#be correct
hdf_dir = os.path.join(script_dir, 'HDF')
csv_dir = os.path.join(script_dir, 'CSV')

#Create the raw save file for the bdot
raw_hdf= hdftools.hdfPath( os.path.join(script_dir, 'run102_PL11B_raw.hdf5') )
lapdtools.readRunProbe(102, 'PL11B', hdf_dir, csv_dir, raw_hdf, verbose=True)

#Create the raw save file for the timing diode
tdiode_hdf = hdftools.hdfPath( os.path.join(script_dir, 'run102_tdiode_raw.hdf5') )
lapdtools.readRunProbe(102, 'tdiode', hdf_dir, csv_dir, tdiode_hdf, verbose=True)

#Combine the two to create a bdot 'full' dataset
full_hdf = hdftools.hdfPath( os.path.join(script_dir, 'run102_PL11B_full.hdf5') )
bdot.bdotRawToFull(raw_hdf, full_hdf, tdiode_hdf=tdiode_hdf, grid=True, verbose=True)
