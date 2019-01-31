# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 12:38:25 2019

@author: Peter
"""

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
     attrs = {}
     
     csvs = getCSVList(csv_dir)
     for csv_file in csvs:
         csvdict = opencsv(csv_file)
         #Experiment runs  
         if csvType(csvdict) == csvType(run=run, probe=None):
             csv_attrs = getRow( csvdict, run=run, probe=None )
             for key in csv_attrs.keys():
                 attrs[key] = csv_attrs[key]
         #Experiment constants
         elif csvType(csvdict) == csvType(run=None, probe=None):
             #experiment constants by definition only have one line b/c
             #they have no row or probe labels
             csv_attrs = csvdict
             for key in csv_attrs.keys():
                 attrs[key] = csv_attrs[key][0]

     #Validate that a matching line was actually found
     #The main program depends on this function to throw an error
     #if this is not the case
     if 'run' in attrs.keys():
             return attrs
     raise ValueError("No spreadsheet found with a row matching input: run=" + str(run))
     
    
    

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
     csvs = getCSVList(csv_dir)
     attrs = {}

     for csv_file in csvs:
         csvdict = opencsv(csv_file)
         if (csvType(csvdict) == csvType(run=run, probe=probe ) 
             or csvType(csvdict) == csvType(run=None, probe=probe )):
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
     