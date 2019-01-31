#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: peter
csvtools.py: Tools for handling metadata CSV files

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



nan_codes = ['NA', '', 'NR', 'NONE']



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
    return csvdict


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





def getRowInd(csvdict, run=None, probe=None):
    
    if not keyExists(csvdict, 'run'):
        run = None
    if not keyExists(csvdict, 'probe'):
        probe = None
    
    #Deal with the case of the experiment constants file
    if csvType(csvdict) == csvType(run=None, probe=None):
        return [0]

    if run is not None and probe is None:
        rowind = [i for i, v in enumerate(csvdict['run']) if
                 v[0] == str(run)]
    elif run is None and probe is not None:
        rowind = [i for i, v in enumerate(csvdict['probe']) if
                 str( v[0] ).lower().strip()  == str(probe).lower().strip()]
    elif run is not None and probe is not None:
        rowind =  [i for i, v in enumerate(csvdict['run']) if 
                 (csvdict['run'][i][0] == str(run) and 
                  str(csvdict['probe'][i][0]).lower().strip()  == str(probe).lower().strip() )  ]
    else:
        rowind = []
    
    #By design there should only be ONE line in a csv file that meets these
    #requirements, so it is safe to conver this to a scalar
    #if the array is empty, the row must not exist.

    if len(rowind) == 0:
        return None
    else:
        return rowind
    
def getRow(csvdict, run=None, probe=None):
    rowind = getRowInd(csvdict, run=run, probe=probe)
    
    if rowind is None:
        return None
    elif len(rowind) > 1:
        raise ValueError("Row specified is not unique! run=" + str(run) 
            + ', probe=' + str(probe) + ' in ' + str(csvdict['csv_filename']))
    else:
        rowind =rowind[0]
        return {key:( fixType(value[rowind][0]) , value[rowind][1] ) 
                for (key,value) in csvdict.items()  }        
    
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
    rowind = getRowInd(csvdict, run=run, probe=probe)
    
    if rowind is None:
        return False
    else:
        return True


def fixType(val):
    """Coerces a value string into an int, str, or float type

    Parameters
    ----------
        val: str

    Returns
    -------
        val: str/ int/ float 
        depending on the inherent type of the variable
    """
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
    return val


def getCSVList(csv_dir):
    """Return a list of all CSV files in a directory. Recursive in file
    structure. Ignores files starting with '.' that are assumed to be hidden.

    Parameters
    ----------
        csv_dir: string 
            directory filepath

    Returns
    -------
        output: list
        List of absolute filepaths to all of the CSV files found.
    """
    output = []
    #Walk the directory tree to find all files
    for root, dir, files in os.walk(csv_dir):
        for f in files:
            name, ext = os.path.splitext(f)
            #If file is a CSV and not hidden (starting with '.')
            if ext == '.csv' and name[0] != '.':
                fpath = os.path.join(root, f)
                output.append(fpath)
    return output






     
def getAllAttrs(csv_dir, run, probe):
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
     csvs = getCSVList(csv_dir)
     attrs = {}

     for csv_file in csvs:
         csvdict = opencsv(csv_file)
         csv_attrs = getRow( csvdict, run=run, probe=probe)
         if csv_attrs is not None:
             for key in csv_attrs.keys():
                 attrs[key] = csv_attrs[key]
             
     #Validate that a matching line was actually found
     #The main program depends on this function to throw an error
     #if this is not the case
     if 'run' in attrs.keys() and 'probe' in attrs.keys():
             return attrs
     raise ValueError("No spreadsheet found with row matching input: run=" 
                      + str(run) +', probe= ' + str(probe))
     
     
     
     
     
    
def getProbeList(csv_dir, run):
    """ Returns a list of all of the probes associated with a given
    run number
    
    Parameters
    ----------
        csv_dir: str
            Path to directory that holds csv files
        run: int
            Run number to find probes for
    Returns
    -------
       list [ (str,str)]
       
       List of tuples (probe name, probe type) associated with this run.
       If no probe_type field is found, defaults to "unknown"
    """
    probelist = []
    probetypes = {}
    
    #Load a list of all the CSV files in the dir
    csvs = getCSVList(csv_dir)
    
    for csv_file in csvs:
        csvdict = opencsv(csv_file)
        
        #If the CSV file is a probe_run file, add all the probes for the run
        #to the probelist
        if  csvType(csvdict) == csvType(run=True, probe=True):
            #Get all the row indices coorespondng to this run
            row = getRowInd(csvdict, run=run)
            if row is not None:
                for r in row:
                    probelist.append(csvdict['probe'][r][0])
                    
                    #If probetypes are defined here in this file, add them to the
                    #dictioanry
                    if keyExists(csvdict, 'probe_type'):
                        probetypes[csvdict['probe'][r][0]] = csvdict['probe_type'][r][0]

        
        #If the CSV file is a probe_constants file, figure out which type each
        #probe is and store it in a dictioanry
        elif csvType(csvdict) == csvType(run=None, probe=True):
            if keyExists(csvdict, 'probe') and keyExists(csvdict, 'probe_type'):
                probes = csvdict['probe']
                for i,p in enumerate(probes):
                    probetypes[p[0]] = csvdict['probe_type'][i][0]
      
    #Convert to a set and back to remove any duplicates         
    probelist = list(set(probelist))
    output = []
    
    #Go through the list and pair probe names with types.
    #If any are not in the dictonary, label them as unknown
    for p in probelist:
        if p in probetypes.keys():
            output.append(  (p, probetypes[p])  )
        else:
            output.append(   (p, 'unknown')   )

    return output



def getRunList(csv_dir):
 #Load a list of all the CSV files in the dir
    csvs = getCSVList(csv_dir)
    
    runlist = []
    for csv_file in csvs:
        csvdict = opencsv(csv_file)
        if keyExists(csvdict, 'run'):
            for entry in csvdict['run']:
                x = fixType(entry[0])
                if x is not np.nan:
                    runlist.append( fixType(entry[0]) )
    
    #Convert to a set and back to remove duplicates
    runlist = list(set(runlist))
    
    
    
def missingKeys(attrs, req_keys, fatal_error=True, missing_keys = []):
    for k in req_keys:
        if not k in attrs.keys():
            missing_keys.append(k)
    if len(missing_keys) > 0:
        if fatal_error:
            raise ValueError("Missing columns in metadata attributes! The following keys were not found, and are required: " + str(missing_keys))
        else:
            return True
    else:
        return False
    
    
    







if __name__ == "__main__":
    fname = r"F:/LAPD_Mar2018/METADATA/CSV/langmuir_runs.csv"
    csv_dir = r"F:/LAPD_Mar2018/METADATA/"

    csvdict = opencsv(fname)
    
    #print( getRow(csvdict, run=102, probe=None) )
    
    #print( findValue( csvdict, 'probe', run=102, probe='PL11B'  ) )
    #v = keyExists(csvdict, 'probe_origin_zz' )
    #print(getAttrs( csvdict,  run=None, probe='PL11B'  ))
    
    #print(rowExists(csvdict, run=50, probe = "LAPD1")  )

    
    #print( findCSV(csv_dir, run=None, probe=None) )
    
    
    #print( getRunLevelAttrs(csv_dir, 80))
    #print( getProbeLevelAttrs(csv_dir,80, 'PL11B'))
    #print( getAllAttrs(csv_dir,80, 'PL11B'))
    
    #print(getProbeList(csv_dir, run=1))

    #print( getRunList(csv_dir) )