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


if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/bdot_runs_LAPD_Mar2018.csv"

    csvobj = opencsv(fname)
    
    print(csvobj.keys())
    print(csvobj['probe'])