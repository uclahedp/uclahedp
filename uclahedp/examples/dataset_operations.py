#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 12:50:55 2019

@author: peter

Example of using dataset routines to process a full datafile
 
"""

from uclahedp.tools import hdf as hdftools
from uclahedp.tools import dataset as dataset

from uclahedp.process import bdot as bdot


#probe = 'LAPD3'
probe = 'LAPD_C2'

src = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "run34_LAPD3_full.hdf5"))
trimmed = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "trimmed.hdf5"))
avgreps = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "avgreps.hdf5"))
thinned = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "run34_" + probe + "_full_thinned.hdf5"))
current = hdftools.hdfPath(os.path.join("/Volumes","PVH_DATA","LAPD_Jul2019","FULL", "run34_" + probe + "_full_current.hdf5"))

print("Trimming time dimension")
dataset.trimDim(src, trimmed, 'time', val_bounds=[0,50e-6], verbose=True)
print("Averaging over reps dimension")
dataset.avgDim(trimmed, avgreps, 'reps', verbose=True, delsrc=True)
print("Thinning time dimension")
dataset.thinBin(avgreps, thinned, 'time', bin=30, verbose=True, delsrc=True)
print("Generate current save file from B")
bdot.fullToCurrent(thinned, current, verbose=False)