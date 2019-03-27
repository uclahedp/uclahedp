#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hrrseq.py: (Non-GUI) Methods to make Terminal "High-Rep-Rate" software sequences

For now, just a test

Created by Scott Feister on Wed Mar 27 06:42:34 2019
"""

import numpy as np
import datetime

class TSequence():
    def __init__(self, nsteps, comment=None):
        self.nsteps = nsteps # Number of steps in this sequence
        self.comment = comment
        self.tdevices = None # List of tdevices

    def write(self, filename):
        """
        Write shot sequence to file readable by HEDP Terminal software.
           Will overwrite if file already exists.
        """
        #TODO: Check integrity of TSequence before trying to write (e.g. all step numbers matching)
        #   - Comments must be strings
        #   - Nsteps must be a positive integer greater than zero
        #   - All sequences must be correct and match with the nsteps
        
        with open(filename, 'w') as f:
            ## Write the header
            if self.comment is not None:
                for line in self.comment.splitlines():
                    f.write("% " + line + "\n")
                f.write("%\n")
            now = datetime.datetime.now()
            f.write("% Sequence file created " + now.strftime("%Y-%m-%d %H:%M") + "\n")
            for tdev in self.tdevs:
                f.write("RESOURCE ID=" + str(tdev.resource_id) + "\t")
            f.write("\n")
            for tdev in self.tdevs:
                f.write("CHANNEL ID=" + str(tdev.channel) + "\t")
            f.write("\n")
            for tdev in self.tdevs:
                if tdev.unit is not None:
                    f.write("UNIT=" + str(tdev.unit).upper() + "\t")
                else:
                    f.write("\t")
            f.write("\n")
            f.write("###\n")
            
            ## Write the values
            for i in range(self.nsteps):
                for tdev in self.tdevs:
                    f.write(str(tdev.vals[i]) + "\t")
                f.write("\n")                

class TDevice():
    """
    One device in the terminal sequence
    """
    def __init__(self, resource_id, channel=0, unit=None):        
        self.unit = unit
        self.resource_id = resource_id
        self.channel = channel
        self.vals = None
    
    def validate(self, nsteps=None):
        """ Validates that the TDevice instance is ready to write to file, etc."""
        if type(self.vals) != np.ndarray:
            raise Exception("Invalid 'vals' must be of type ndarray")
                
    
        
if __name__ == "__main__":
    # Initialize sequence
    nsteps = 50
    nreps = 5
    comment = "First test of lens stage in an HRR sequence"
    comment += "\nAuthor: Scott Feister"
    
    tseq = TSequence(nsteps*nreps, comment=comment)

    # Define devices and step values
    lens_stage = TDevice(146, unit="mm", channel=0)
    lens_stage.vals = np.linspace(10, 70, nsteps) # Sequence of values
    lens_stage.vals = np.tile(lens_stage.vals, [nreps,1]).T.flatten() # Repeated sequence of values
        
    # Place devices into the sequence
    tseq.tdevs = [lens_stage]
    
    # Write the sequence to file
    tseq.write(r"C:\Users\scott\Documents\temp\march2019\test.txt")
    
    

        
    
