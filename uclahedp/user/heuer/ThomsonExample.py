#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:22:23 2020

@author: peter
"""

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit as curve_fit

from uclahedp.analyze.Thomson import ThomsonScattering as TS



def inst_fcn(wavelength):
    #Create an instrument function that represents the spectral response
    #of the spectrometer recording the Thomson spectrum
    res= 0.1
    fcn = np.exp(-0.5*np.power(wavelength/res,2))/(res*np.sqrt(2*np.pi))
    return fcn



fig, ax = plt.subplots()


#***********************************
#Setup parameters for the synthetic dataset.
#***********************************
Te  = 5e-3 #keV
fract=np.array([1.0])
Ti = 1e-3 #keV
Z = np.array([1])
Mu = np.array([4])
ne = 1e13 #cm^-3

laser_wavelength = 532 #3.5e15 rad/s 532 nm
wvlgth_range=(500,600,.01)

#Angles
sa = 63
dphi =90
     
vpar=1e6
vperp=0
eidrift=0 
gamma=0 #Angled (deg) between k and edrift

blockw = 0.4

#Create the wavelength array
wavelength = np.arange(wvlgth_range[0], wvlgth_range[1], wvlgth_range[2])

#Create some synthetic data, including using the provided insturment fcn
synthdata = TS.formFactor(Te, fract, Ti, Z, Mu, ne,
                   laser_wavelength, sa, dphi,
                   vpar=vpar, vperp=vperp, eidrift=eidrift, gamma=gamma,
                   wvlgth_range=wvlgth_range,
                   blockw=None, inst_fcn=inst_fcn)

#Add synthetic noise on top of the data
synthdata += (np.random.random(synthdata.size)-0.5)*3e-23

xrng = 10
ax.set_xlim(laser_wavelength-xrng, laser_wavelength+xrng)
ax.plot(wavelength, synthdata)






#***********************************
#Fitting Routine
#***********************************

#Block the ion feature (sharp central peak) in the data
#This is important for fitting the electron feature, because otherwise the ion
#feature will dominate the fit. The width chosen depends on the width of hte ion feature
#Obviously it's crucial that the region blocked is the same in the data
#and the fit fcn so it doesn't just fit this dip...
synthdata =  np.where(np.abs(wavelength - laser_wavelength)< blockw, 0, synthdata)


#Setup a curve fitting function as a lambda function of the form factor fcn.
#we'll choose the electron temperature and density as the fit parameters
#The rest of the parameters we will assume are known, and we'll set their values
#explicitly. We'll use the same wavelength range and insturment fcn as the synth
#data. 

#We'll also block the probe laser feature so the fit is dominated by the wings
#Blocking the central ion feature is very important to allowing the fit to
#converge on the electron wings

fitfcn = lambda wavelength, TeVAR, neVAR : TS.formFactor(TeVAR, fract, Ti, Z, Mu, neVAR,
                   laser_wavelength, sa, dphi,
                   vpar=vpar, vperp=vperp, eidrift=eidrift, gamma=gamma,
                   wvlgth_range=wvlgth_range,
                   blockw=blockw, inst_fcn=inst_fcn)




#Set the initial guess values for the fit parameters
#Intentionally making these slightly wrong
p0 = np.array([300e-3, 8e13])

popt, pconv = curve_fit(fitfcn, wavelength, synthdata, p0=p0)

print(popt)

#Extract the estimated values from the fit
TeFit, neFit = popt

#Make a a curve from the fit function at the best-fit parameters
fitcurve = fitfcn(wavelength, TeFit, neFit)

#Plot the best-fit curve
ax.plot(wavelength, fitcurve)






