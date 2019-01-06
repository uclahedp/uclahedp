#Example analyisis 

import os
import uclahedp.bapsftools as bt
import uclahedp.hedpConstants as c
from uclahedp.sdfClass import sdfarr, sdfattrs
import astropy.units as u

exp_dir = r'C:\Users\Peter\Desktop\example01\EXPERIMENT_DIR'

sname = r"C:\Users\Peter\Desktop\testsave2.h5"


run = 56
probe = 'LAPD1'

# This is the absolute path to the bdot_runs_csv, using constants stored in the package
bdot_runs_csv= exp_dir + c.bdot_runs_csv

# Read the data in from the HDF file using the bdot_runs_csv for information
obj = bt.readRunProbe(56, 'LAPD1',  exp_dir, bdot_runs_csv)

# Save the file 
savefile = exp_dir  + c.bdot_dir  +  'run' + str(run) + '_' + str(probe) + '.h5'
obj.saveHDF(savefile, hdfpath = '/Run56/LAPD1/')

#Plot the raw signal
print(obj.labels)
print(obj.data.shape)
obj.collapseDim('channels', 1)
obj.avgDim('reps')
obj.collapseDim('x', 0)
obj.convertAxisUnit('t', u.us)
obj.plot(xrange=[0,20])


