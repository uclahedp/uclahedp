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
    """ Opens a CSV files as a dictionary
    Keys are taken from the first row of the csv. Columns without
    a key are ignored
    Skips the 2nd and 3rd lines of the CSV, which are human-readable titles

    Parameters
    ----------
        fname : str
            Filepath of the CSV file

    Returns
    -------
        csvdict: dict
            dictionary of key:list pairs, where each key is a header from
            row 1 and each list contains entries from rows [4,-1]
    """
    csvdict = defaultdict(lambda: [])
    with open(fname, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        keys = reader.fieldnames
        # Skip the first two lines (human-readable comments)
        next(reader)
        next(reader)
        for row in reader:
            for key in keys:
                csvdict[key].append(row[key])
    return csvdict


def findvalue(csvdict, key, run=None, probe=None):
    """ Picks out a value from a csv file dictionary
    Returns a list of all values which match the given key, run, and probe.

    Parameters
    ----------
        csvdict: dict
            Dictonary object to search
        key: str
            Key for the value you are searching for
        run: int
            Integer run number for the value you are searching for
        probe: str
            Name of the probe you are searching for.

    Returns
    -------
        value: list
    """
    if probe is None:
        value = [v for i, v in enumerate(csvdict[key]) if
                 csvdict['run'][i] == str(run)]
    elif run is None:
        value = [v for i, v in enumerate(csvdict[key]) if
                 csvdict['probe'][i] == str(probe)]
    else:
        value = [v for i, v in enumerate(csvdict[key]) if
                 (csvdict['run'][i] == str(run) and
                  csvdict['probe'][i] == str(probe))]
    
    if len(value) == 0:
        return []
    elif len(value) == 1:
        value = value[0]
    
    
    try: 
        value = float(value)
        
        if value.is_integer():
            return int(value)
        else:
            return value
        
    except ValueError:
        return str(value)




if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/bdot_runs.csv"

    csvdict = opencsv(fname)
    v = findvalue(csvdict, 'probe_origin_z', run=40, probe='LAPD4')
    print(v)
    print(type(v))
