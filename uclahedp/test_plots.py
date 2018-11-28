import matplotlib.pyplot as plt
import numpy as np
import h5py
import rawtools

fname_h5 = r"/Volumes/PVH_DATA/LAPD_Mar2018/RAW/run56_LAPD1_pos_raw.h5"

with h5py.File(fname_h5, 'r') as f:
    gridded, xgv, ygv, zgv = rawtools.regrid(f)
    time = 1e6*rawtools.timevector(f, pos=35, rep=0)
    dbtrace = gridded[:, 35, 0, 0, 0, 1]

btrace = np.cumsum(dbtrace)

plt.plot(time, btrace)
plt.axis([0, 40, -3, 3])
plt.xlabel('t (us)')
plt.ylabel('BY (arb)')
plt.show()
