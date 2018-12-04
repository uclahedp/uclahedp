#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hedpConstants.py: This file holds a list of string constants
These are used throughout the analysis code to avoid hardcoding in file names
@author: peter
"""
#
version = '1.0.0'
#FILE PATHS
#All filepaths and directories are given from a base directory
raw_dir = '/RAW/'
hdf_dir = '/HDF/'
metadata_dir = '/METADATA/CSV/'

tdiode_dir = '/TDIODE/'
bdot_dir = '/BDOT/'
langmuir_dir = '/LANGMUIR/'
interferometer_dir = '/INTERFEROMETER/'

#FILES
main_runs_csv = metadata_dir + 'main_runs.csv'

tdiode_runs_csv = metadata_dir +  'tdiode_runs.csv'
bdot_runs_csv = metadata_dir +  'bdot_runs.csv'
bdot_probes_csv = metadata_dir +  'bdot_probes.csv'
langmuir_runs_csv = metadata_dir +  'langmuir_runs.csv'
langmuir_probes_csv = metadata_dir +  'langmuir_probes.csv'
interferometer_runs_csv = metadata_dir + 'interferometer_runs.csv'
