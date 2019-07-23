#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hrrseq.py: (Non-GUI) Methods to make Terminal "High-Rep-Rate" software sequences

Created by Scott Feister on Wed Mar 27 06:42:34 2019
"""

import numpy as np
import datetime

class TSequence():
    """
    A terminal sequence.
        comment: Anything you'd like to note regarding this sequence.
        nsteps: Number of steps in the sequence, excluding repetitions below
        nreps:  Number of repetitions of DAQ at each step (e.g. to build statistics),
                    implemented here by duplicating; total acquisitions is nsteps x nreps.
        tdevs:  A list of objects of class TDevice.
    """
    
    def __init__(self, nsteps, nreps=1, comment=None):
        self.nsteps = nsteps # Number of steps in this sequence
        self.comment = comment
        self.nreps = nreps
        self.tdevs = None # List of tdevices

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
            f.write("% Steps: " + str(self.nsteps) + ", Reps/step: " + str(self.nreps) + ".\n")
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
            
            ## Write the values, repeating nreps times at each step
            for i in range(self.nsteps):
                for j in range(self.nreps):
                    for tdev in self.tdevs:
                        f.write(str(tdev.vals[i]) + "\t")
                    f.write("\n")                

    def validate(self):
        # TODO: Make this section more failproof, and errors more clear on next steps
        if self.tdevs is None:
            raise Exception("Failed validation. No devices defined.")
        for tdev in self.tdevs:
            tdev.validate(nsteps=self.nsteps)
            
class TDevice():
    """
    One device in the terminal sequence
    """
    def __init__(self, resource_id, channel=0, unit=None):        
        self.unit = unit
        self.resource_id = resource_id
        self.channel = channel
        self.vals = None
    
    def __str__(self):
        return str("Resource ID " + str(self.resource_id) + ", Channel " + str(self.channel))
    
    def validate(self, nsteps=None):
        """ Validates that the TDevice instance is ready to write to file, etc."""
        if type(self.vals) != np.ndarray:
            raise Exception("Failed validation. 'vals' must be of type ndarray for device " + str(self))
        if nsteps is not None and len(self.vals) != nsteps:
            raise Exception("Failed validation. Number of steps mismatch for device " + str(self))
        
    
        
if __name__ == "__main__":
    # Initialize sequence
    nsteps = 21
    nreps = 1
    comment = "Coarse Scan of LAPD Lens, Peter's July 2019 Experiment"
    comment += "\nAuthor: Scott Feister"
    tseq = TSequence(nsteps, nreps=nreps, comment=comment)

    # Define devices and step values
    lens = TDevice(145, unit="mm", channel=0)    
    lens.vals = np.linspace(260, 660, tseq.nsteps) # Sequence of values

    # Safety check for LAPD lens limits of 0 mm, 740 mm
    if np.max(lens.vals) > 740 or np.min(lens.vals) < 0:
        raise Exception("Outside known acceptable range of LAPD Lens.")
        
    # Place devices into the sequence
    tseq.tdevs = [lens]
    
    # Write the sequence to file
    tseq.write(r"C:\Users\scott\Documents\temp\jul2019\test.txt")
    
    

        
    
