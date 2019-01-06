#Example analyisis 

import os
import uclahedp.bapsftools as bt
import uclahedp.hedpConstants as c

exp_dir = r'C:\Users\Peter\Desktop\example01\EXPERIMENT_DIR'


run = 56
probe = 'LAPD1'



bdot_runs_csv= exp_dir + c.bdot_runs_csv


obj = bt.readRunProbe(56, 'LAPD1',  exp_dir, bdot_runs_csv)

print(exp_dir)
savefile = exp_dir  + c.bdot_dir  +  'run' + str(run) + '_' + str(probe) + '.h5'
print(savefile)

obj.saveHDF(savefile, hdfpath = 'test')