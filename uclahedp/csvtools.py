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


def findvalue(csvdict, key, run=0, probe=None):
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
    else:
        value = [v for i, v in enumerate(csvdict[key]) if
                 (csvdict['run'][i] == str(run) and
                  csvdict['probe'][i] == str(probe))]
    return value


if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/bdot_runs_LAPD_Mar2018.csv"

    csvdict = opencsv(fname)
    print(findvalue(csvdict, 'probe_xpos', run=40, probe='LAPD1'))
