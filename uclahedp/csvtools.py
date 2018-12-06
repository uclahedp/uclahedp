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


def keyExists(csvdict, key):
    """ Tests to see if a key exists in a csv dictionary
    Parameters
    ----------
        key: str
            Key for the value you are searching for
            
        csvdict: dict
            Dictonary object to search

    Returns
    -------
        check: boolean
    """
    return key.lower().strip() in [(k.lower().strip()) for k in csvdict.keys()]



def getAttrs(csvdict, run=None, probe=None):
    """ Returns a dictionary of all of the keyed csv lines for a run/probe/both
    
    Parameters
    ----------
        csvdict: dict
            Dictonary object to search
        run: int
            Integer run number for the value you are searching for
        probe: str
            Name of the probe you are searching for.

    Returns
    -------
        attrs: dict
    """
    attrs = {}
    for key in csvdict.keys():
        attrs[key] = findValue(csvdict, key, run=run, probe=probe)
    
    return attrs
    


def findValue(csvdict, key, run=None, probe=None):
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
    
    # These values all code to "None" 
    int_none = [-111, -999]
    float_none = [float(i) for i in int_none]
    str_none = ['none', 'na']
    
    if probe is None:
        value = [v for i, v in enumerate(csvdict[key]) if
                 csvdict['run'][i] == str(run)]
    elif run is None:
        value = [v for i, v in enumerate(csvdict[key]) if
                 str( csvdict['probe'][i] ).lower().strip()  == str(probe).lower().strip()]
    else:
        value = [v for i, v in enumerate(csvdict[key]) if 
                 (csvdict['run'][i] == str(run) and 
                  str(csvdict['probe'][i]).lower().strip()  == str(probe).lower().strip() )  ]
    
    if len(value) == 0:
        return None
    elif len(value) == 1:
        value = value[0]
    
    
    try: 
        value = float(value)
        
        #IS AN INT
        if value.is_integer():
            if int(value) in int_none:
                return None
            else:
                return int(value)
        
        #IS A FLOAT
        else:
            if float(value) in float_none:
                return None
            else:
                return value
   
    except ValueError:
        #IS A STRING
        if str(value).lower() in str_none:
            return None
        else:
            return str(value).strip()




if __name__ == "__main__":
    fname = r"/Volumes/PVH_DATA/LAPD_Mar2018/METADATA/CSV/bdot_runs.csv"

    csvdict = opencsv(fname)
    v = findValue( csvdict, 'probe_origin_z', run=50, probe='lapd1'  )
    print(v)
    v = keyExists(csvdict, 'probe_origin_zz' )
    print(v)
    attrs = getAttrs( csvdict,  run=50, probe='lapd1'  )
    print(attrs)