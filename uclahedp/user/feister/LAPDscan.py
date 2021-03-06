#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simplescan.py: Simple 1D translation sequence for the LAPD Lens

For Peter's July 2019 Experiment

Created by Scott Feister on Tue Jul 23 14:15:40 2019
"""

import numpy as np
import datetime
import os
from hrrseq import TSequence, TDevice

if __name__ == "__main__":
    # Initialize sequence
    
    ds = 10
    dstart = 0
    dstop = 720
    
    npos = int((dstop-dstart)/ds)+1
  
    nreps = 10
    comment = "Coarse Scan of LAPD Lens, Peter's July 2019 Experiment"
    comment += "\nAuthor: Scott Feister"
    tseq = TSequence(npos, nreps=1, comment=comment)


    lens = TDevice(145, unit="mm", channel=0)    
    #lens.vals = np.linspace(260, 660, tseq.nsteps) # Sequence of values
    lens.vals = np.arange(dstart, dstop+ds, ds)
    
    
    npos = lens.vals.size
    print("Positions: " + str(npos))
    print("Minutes: " + str(npos*nreps/60))
        
    # Safety check for LAPD lens limits of 0 mm, 740 mm
    if np.max(lens.vals) > 740 or np.min(lens.vals) < 0:
        raise Exception("Outside known acceptable range of LAPD Lens.")
        
    # Place devices into the sequence
    tseq.tdevs = [lens]
    
    # Write the sequence to file
    outdir = r"PVH_LAPD_Outputs"
    datestr = datetime.datetime.now().strftime("%Y_%m_%d")
    outname = datestr + "_LAPD_Lens_" + str(dstart) + "_to_" + str(dstop) \
                + "_in_" + str(npos) + "_steps_of_" + str(nreps) + "_reps.txt"
    outpath = os.path.join(outdir, outname)
    tseq.write(outpath)
    
    

        
    
