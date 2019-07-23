#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simplescan.py: Simple 1D translation sequence for the LAPD Lens

For Peter's July 2019 Experiment

Created by Scott Feister on Tue Jul 23 14:15:40 2019
"""

import numpy as np
from hrrseq import TSequence, TDevice

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
    
    

        
    
