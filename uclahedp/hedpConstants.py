#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hedpConstants.py: This file holds a list of string constants
These are used throughout the analysis code to avoid hardcoding
@author: peter
"""
from os.path import join
#
version = '1.0.0'
#FILE PATHS
#All filepaths and directories are given from a base directory
raw_dir = '/RAW/'
hdf_dir = '/HDF/'
metadata_dir = '/METADATA/'

tdiode_dir = '/TDIODE/'
bdot_dir = '/BDOT/'
langmuir_dir = '/LANGMUIR/'
interferometer_dir = '/INTERFEROMETER/'

not_applicable_codes = ['NA']
not_recorded_codes = ['', 'NR']