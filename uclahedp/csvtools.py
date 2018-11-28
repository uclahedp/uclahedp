#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
csvtools.py: Tools for handling metadata CSV files

Created on Wed Nov 28 12:01:01 2018

@author: peter
"""

import csv
from collections import defaultdict


def opencsv(fname):
    d = defaultdict(lambda: [])
    with open(fname, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        keys = reader.fieldnames
        # Skip the first two lines (human-readable comments)
        next(reader)
        next(reader)
        for row in reader:
            for key in keys:
                d[key].append(row[key])
    return d


def findvalue(csvdict, key, run=0, probe=None):

    if probe == None:
        value = [v for i,v in enumerate(csvdict[key]) if 
                 csvdict['run'][i] == str(run)]
    else:
        value = [v for i,v in enumerate(csvdict[key]) if 
                 (csvdict['run'][i] == str(run) and 
                  csvdict['probe'][i] == str(probe)  )]
    return value


if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/bdot_runs_LAPD_Mar2018.csv"

    csvdict = opencsv(fname)
    print(findvalue(csvdict, 'probe_xpos'))