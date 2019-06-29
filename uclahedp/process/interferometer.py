# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 11:33:00 2019

@author: Peter
"""

import os, h5py
from uclahedp.tools import hdf as hdftools
from uclahedp.tools import dataset

import numpy as np
import matplotlib.pyplot as plt


def abs_phase(sig, ref):
     #Calculate the ffts
     ref_fft = np.fft.fft(ref)
     sig_fft = np.fft.fft(sig)
     nfreq = ref_fft.shape[0]
     
     #Determine the dominant frequency from the reference signal
     nyquist = int(nfreq/2)
     ref_freq = np.argmax(ref_fft[1:nyquist])

     #Apply bandpass
     #Intentionally zeros negative frequencies 
     nfreq = ref.shape[0]
     delta = int(nfreq*0.01) #Arbitrary small window of frequencies
     ref_fft[0:ref_freq-delta] = 0
     ref_fft[ref_freq+delta:-1] = 0
     sig_fft[0:ref_freq-delta] = 0
     sig_fft[ref_freq+delta:-1] = 0

     #Reverse FFT
     sig = np.fft.ifft(sig_fft)
     ref = np.fft.ifft(ref_fft)
     
     #Compute the phase dufference by multiplying by the conjugate
     sig = sig*np.conj(ref)
     #Separate out the phase
     phase = np.angle(sig)
     #Unwrap the phase (correct jumps of more than pi)
     phase = np.unwrap(phase)

     return phase

     
     
     
     
     





if __name__ == "__main__":
     
     sig_file = hdftools.hdfPath(os.path.join("F:","LAPD_Mar2018","RAW", "run83_interferometer_signal.hdf5"))
     ref_file = hdftools.hdfPath(os.path.join("F:","LAPD_Mar2018","RAW", "run83_interferometer_reference.hdf5"))
     dest_file = hdftools.hdfPath(os.path.join("F:","LAPD_Mar2018","FULL", "run83_interferometer.hdf5"))
     
     with h5py.File(sig_file.file) as f:
          sig = f['data'][:,0]
          time = f['time'][:]
          nti = time.shape[0]
          
     with h5py.File(ref_file.file) as f:
          ref = f['data'][0:nti,0]
     
     phase = abs_phase(sig, ref)
     
     axes = [{'ax':time, 'name':'time', 'unit':'s'}]
     dataset.createDataset(phase/3.14, axes, dest_file, dataunit='rad', attrs=None)

