#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
csvtools.py: Tools for handling metadata CSV files

Created on Wed Nov 28 12:01:01 2018

@author: peter

**CSV types**
Four types of CSV files are supported, differentiated by the presence of 'run'
and 'probe columns.

Experiment constants -> No runs column, no probe column. Data that applies to
every run of the experiment (at the run level)

Experiment runs -> Runs column, no probe column. Data about each run that
applies to every probe  used. Ex. background B-field, gas, etc.

Probe constants -> No runs column, probes column. Data about each probe that is
constant throughout the experiment. Ex. tip area.

Probe runs -> Runs column, probe column. Data about each probe for each run.
Can contain multiple probe rows for each run (but each run,probe pair should be
an unique row.)



"""

import csv
from collections import defaultdict
import os
import numpy as np



nan_codes = ['NA', '', 'NR']



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
        unitdict = dict( next(reader) )#save the units as a dictionary
        next(reader) #Skip the human-readable titles line
        for row in reader:
            for key in keys:
                #Don't include lines with no header key
                if key != '':
                    csvdict[key].append( (row[key], unitdict[key] )    )  
        #csvdict['units'] = unitdict
    return dict(csvdict)


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




def rowExists(csvdict, run=None, probe=None):
    """Checks to see if the given csv file has columns that match a run/probe
    pair

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
        bool: A matching row exists, or it does not.
    """
    if run and not probe :
        val = findValue(csvdict, 'run', run=run)
    elif not run and probe:
        val = findValue(csvdict, 'probe', probe=probe)
    elif run and probe:
        val = findValue(csvdict, 'run', run=run, probe=probe)
    else:
        #Without a run or probe column, this csv file cannot be used as valid
        #metadata
        return False
    
    return bool(val)




def csvType(*args, run=None, probe=None):
    """ Determines the type of a csvdict in one of two ways
    
    If an arg is supplied, it is assumed to be a csvdict. The function will 
    check for run and probe keys, and report the type of that csv dict
    
    
    Alternately, if the run and probe keywords are set, the program will return
    what type of csv file WOULD be compatible with that combination of keys.

    Parameters
    ----------
        csvdict: dict OPTIONAL
            Dictonary object to search
        run: int/None or Boolean
        
        probe: str/None or Boolean
    Returns
    -------
        str:
            Type of csv file
    """
    
    #If an argument is supplied, assume it is a csv dict
    if len(args) == 1:
        csvdict = args[0]
        run = keyExists(csvdict, 'run')
        probe = keyExists(csvdict, 'probe')

    if not run and not probe:
        return 'experiment_constants'
    elif run and not probe:
        return 'experiment_runs'
    elif not run and probe:
        return 'probe_constants'
    elif run and probe:
        return 'probe_runs'
    
    
    


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
    
    #If the csvdict doesn't have the right key, it cannot have the data!
    if run and not keyExists(csvdict, 'run'):
        return None
    if probe and not keyExists(csvdict, 'probe'):
        return None

    if run is None and probe is None:
         arr = [v for v in csvdict[key]]
    elif run is not None and probe is None:
        arr  = [v for i, v in enumerate(csvdict[key]) if
                 csvdict['run'][i][0] == str(run)]
    elif run is None and probe is not None:
        arr  = [v for i, v in enumerate(csvdict[key]) if
                 str( csvdict['probe'][i][0] ).lower().strip()  == str(probe).lower().strip()]
    elif run is not None and probe is not None:
        arr  = [v for i, v in enumerate(csvdict[key]) if 
                 (csvdict['run'][i][0] == str(run) and 
                  str(csvdict['probe'][i][0]).lower().strip()  == str(probe).lower().strip() )  ]
                 

    if len(arr) == 0:
        return None
    elif len(arr ) >= 1:
        arr  = arr[0]
    
    #Decompose the tuple into a value and a unit string
    val = arr[0]
    unit = str(arr[1])
        
    #Determine if data is an int, float, or string, and type it accordingly
    try: 
        val = float(val)
        #IS AN INT
        if val.is_integer():
                val= int(val)
    except ValueError:
        #IS A STRING
        if val in nan_codes:
            val = np.NAN
        else:
            val = str(val).strip()
    
    return (val, unit)
            
        
        
        
        
        
        
        
        
def findCSV(csv_dir, run=None, probe=None):
    """ Find CSV files in a directory that coorespond to a run/probe pair
    The run/probe pair is exclusionary: searching for just a run will only
    return csv files with JUST a run column, not those with a run and probe
    column
    
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
       list: List of absolute file paths to csv files matching request
    """
    
    #Determine the type of CSV file implied by the keword arguments
    ctype = csvType(run=run, probe=probe)

    output = []
    
    #Walk the directory tree to find all files
    for root, dir, files in os.walk(csv_dir):
        for f in files:
            name, ext = os.path.splitext(f)
            #If file is a CSV and not hidden (starting with '.')
            if ext == '.csv' and name[0] != '.':
                fpath = os.path.join(root, f)
                csvdict = opencsv(fpath)
                #If the CSV is of the type implied by the run and probe keywords
                if csvType(csvdict) == ctype:  
                    #If the CSV file contains the run/probe requested
                    #Second part of conditional deals with the case of experimental constants data
                    if rowExists(csvdict, run=run, probe=probe) or (run is None and probe is None):
                        output.append(fpath)
    return output
 



def getRunLevelAttrs(csv_dir, run):
     """ Returns a dictionary of all attributes from a directory of CSV files
     that coorespond to a particular run, at the run level. These include all
     attributes about the entire experiment, but not about the probes
    
    Parameters
    ----------
        csv_dir: str
            Path to directory that holds csv files
        run: int
            Name of the probe you are searching for.
    Returns
    -------
       dict: combined attribute dictionary describing the run
       """
     csvs = findCSV(csv_dir, run=run, probe=None) + findCSV(csv_dir, run=None, probe=None)
     attrs = {}
     for csv_file in csvs:
        print(csv_file)
        csv_attrs = getAttrs( opencsv(csv_file ), run=run, probe=None )
        for key in csv_attrs.keys():
            attrs[key] = csv_attrs[key]
     return attrs
    
    

def getProbeLevelAttrs(csv_dir, run, probe):
     """ Returns a dictionary of all attributes from a directory of CSV files
     that coorespond to a particular probe and run. These include attributes
     that describe the probe for that run, as well as constants about the
     probe.
    
    Parameters
    ----------
        csv_dir: str
            Path to directory that holds csv files
        probe: str
            Integer run number for the value you are searching for
    Returns
    -------
       dict: combined attribute dictionary describing the run and probe
    """
     csvs = findCSV(csv_dir, run=None, probe=probe) + findCSV(csv_dir, run=run, probe=probe)
     attrs = {}

     for csv_file in csvs:
         csv_attrs = getAttrs( opencsv(csv_file ), run=run, probe=probe)
         for key in csv_attrs.keys():
             attrs[key] = csv_attrs[key]
     return attrs 
    
    


if __name__ == "__main__":
    fname = r"F:/LAPD_Mar2018/METADATA/CSV/bdot_runs.csv"
    csv_dir = r"F:/LAPD_Mar2018/METADATA/"

    csvdict = opencsv(fname)
    
    #print( findValue( csvdict, 'probe', run=102, probe='PL11B'  ) )
    #v = keyExists(csvdict, 'probe_origin_zz' )
    #print(getAttrs( csvdict,  run=102, probe='PL11B'  ))
    
    #print(rowExists(csvdict, run=50, probe = "LAPD1")  )

    
    #print( findCSV(csv_dir, run=None, probe=None) )
    
    
    print( getRunLevelAttrs(csv_dir, 102))
    #print( getProbeLevelAttrs(csv_dir,102, 'PL11B'))
    