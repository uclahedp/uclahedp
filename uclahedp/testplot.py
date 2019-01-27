# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import h5py

import hdftools
import os

src = hdftools.hdfPath( os.path.join("F:", "LAPD_Mar2018", "RAW", "test_PL11B_full.hdf5") )

tind =600

with h5py.File(src.file, "r") as f:
    print(f[src.group]['data'].shape)
    data = f[src.group]['data'][tind,:,0,:,0,2] + 300
    v = f[src.group]['xaxes'][:]
    h = f[src.group]['zaxes'][:]
    print(data.shape)
    print(h.shape)
    print(v.shape)

    print(f[src.group]['time'][tind]*1e6)
    
    cmin = 0
    cmax = np.max(data)
    nlevels = 100
    levels = np.linspace(cmin, cmax, num=nlevels)


    fig, ax = plt.subplots()
    c= ax.contourf(h, v, data, levels, cmap=plt.cm.RdBu_r)

    cbar = plt.colorbar(c)
    plt.show()
