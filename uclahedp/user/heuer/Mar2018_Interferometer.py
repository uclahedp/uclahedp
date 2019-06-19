# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 17:31:23 2019
This program processes interferometer data for the Mar2018 run
"""
import numpy as np
import os

from uclahedp.tools import csv as csvtools
from uclahedp.tools import hdf as hdftools
from uclahedp.load import ascii


def make_raw():
     probe = 'TonyInterferometer'
     
     
     data_dir = os.path.join("F:","LAPD_Mar2018","INTERFEROMETER")
     out_dir = os.path.join("F:","LAPD_Mar2018","RAW")
     csv_dir = os.path.join("F:","LAPD_Mar2018","METADATA", "CSV")
     metadata_csv = os.path.join(csv_dir, "interferometer_runs.csv")
     
     attrs = csvtools.opencsv(metadata_csv )
     
     runs = attrs['run']
     ref_file = attrs['ref_file']
     sig_file = attrs['sig_file']
     index = attrs['index']
     
     
     for i, r in enumerate(runs):
          if index[i][0] == '0':
               run = int(runs[i][0])
               print("Processing run " + str(run))
               sf = hdftools.hdfPath(os.path.join(data_dir, sig_file[i][0] + '.txt'))
               filename = 'run' + str(run) + '_interferometer_signal.hdf5'
               dest = hdftools.hdfPath(os.path.join(out_dir, filename))
           
          
               ascii.asciiToRaw(sf, dest, skip_header=5, ax=0, axis_name='time', 
                          axis_unit='s', data_unit='V',
                          csv_dir=csv_dir, run=run, probe=probe)
               
               
               if ref_file[i][0] != 'NA':
                    sf = hdftools.hdfPath(os.path.join(data_dir, ref_file[i][0] + '.txt'))
                    filename = 'run' + str(runs[i][0]) + '_interferometer_reference.hdf5'
                    dest = hdftools.hdfPath(os.path.join(out_dir, filename))
          
               ascii.asciiToRaw(sf, dest, skip_header=5, ax=0, axis_name='time', 
                          axis_unit='s', data_unit='V',
                          csv_dir=csv_dir, run=run, probe=probe)



if __name__ == '__main__':
     
     make_raw()