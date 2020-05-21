# -*- coding: utf-8 -*-
"""
thomson_scattering.py 

@author: Peter Heuer

This module contains code for the analysis of Thomson scattering spectra
from plasmas with multiple ion species.

This function was adapted from a program provided in the book
     "Plasma Scattering of Electromagnetic Radiation" by Froula et al.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func_deriv as ZPrime




def spectral_density(wavelength, Te=None, Ti=None, ion_fract=np.array([1.0]), ion_Z=None, 
               ion_Mu=None, ne=None, probe_wavelength=None,
               sa = None, vpar=0.0, vperp=0.0, eidrift=0.0, gamma=0,
               vol=None, solid_angle = 4*np.pi) :
     
     """
     Te (float): Electron temperature in KeV
     
     ion_fract (float ndarray): Relative fractions of each ion species
     present. Must sum to 1.0.
     
     Ti (float): Ion temperature in KeV (all species must have the same temp?)
     
     Z (float ndarray): Ion charge states, normalized to the proton fundamental
     charge
     
     Mu (float ndarray): Ion atomic masses, normalized to the proton mass
     
     ne (float): Plasma electron number density in cm^-3. 
     
     laser_wavelength (float): Incident laser wavelength in nm
     
     sa (float): Scattering angle in degrees
     
     dphi (float): Angle between the polariazation-plane of the incident laser
     and the scattering plane in degrees.
     
     vpar (float): Plasma fluid velocity perpendicular to the probe laser
     in cm/s
     
     vperp (float): Plasma fluid velocity parallel to the probe laser in
     cm/s 
     
     eidrift (float ndarray): Relative drift between the electrons and each ion
     species in cm/s
     
     gamma (float): Angle (in degrees) between k and the drift velocity eidrift
     
     
     blockw (None or float): If not none, block a range of width blockw
     around the laser wavelength to zero
     
     inst_fcn (None or function): If not none, represents an insturment
     function (as a python function that takes a wavelength argument and
     returns a normalized insturment function)
     
     
     """
     
     #TODO: Support multiple ion populations with different vpar and vperp
     
     if np.sum(ion_fract) != 1.0:
          print("WARNING: Sum(IonFractions) != 1.0. Normalizing array." )
          ion_fract = ion_fract/np.sum(ion_fract)
          
    
     #Define Constants
     C=2.99792458e10 #Speed of light, cm/s
     re=2.8179e-13 #Classical Electron radius in cm
     Me=5.6e-19 # Electron mass in keV / (cm/s)^2
     Mp=Me*1836.1 #Proton mass
     Mi=ion_Mu*Mp #Mass of each ion species
     

     #Calculate some plasma values
     Esq = Me*C**2*re
     const = np.sqrt(4*np.pi*Esq/Me)
     
     mean_Z = np.sum(ion_Z*ion_fract) #dimensionless
     ni = ion_fract*ne/mean_Z #cm^-3
     wpe = const*np.sqrt(ne) #Electron plasma frequency in rad/s
     wpi = const*ion_Z*np.sqrt(ni*Me/Mi) # rad/s
     vTe = np.sqrt(Te/Me) #cm/s
     vTi = np.sqrt(Ti/Mi) #cm/s
    
     #Convert all the angles to radians
     sa = sa*2*np.pi/360.0 #rad
     #dphi = dphi*2*np.pi/360.0 #rad
     gamma = gamma*2*np.pi/360.0 #rad
     
     #Convert wavelengths to angular frequencies (electromagnetic waves, so 
     #phase speed is c)
     ws = 2*np.pi*C*1e7/wavelength #rad/s
     wl = 2*np.pi*C*1e7/probe_wavelength #rad/s

     #w is the frequency SHIFT
     w = ws - wl #Eq. 1.7.8 in Sheffield
     
    
     #Wavenumbers in the plasma (wpe emerges from plasma index of refraction?)
     ks = np.sqrt(ws**2 - wpe**2)/C #Eq. 5.4.1 in Sheffield, units 1/cm
     kl = np.sqrt(wl**2 - wpe**2)/C #Eq. 5.4.2 in Sheffield, units 1/cm
     
     #k is the frequency SHIFT (as illustrated in Fig.1.5 of Sheffield)
     k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(sa)) #Eq. 1.7.10 in Sheffield, units 1/cm
     
     
     #Compute the wave frequency Doppler shifted into the plasma rest frame
     kdotv = (kl - ks*np.cos(sa))*vpar - ks*np.sin(sa)*vperp # Units rad/s
     wdoppler = w - kdotv #Units rad/s
     #TODO change this to compute a different wdoppler for each ion species, 
     #to allow for each ion population to have its own velocity?
     
     
     #These dimensoonless constants correspond to alpha in Schaeffer's thesis, 
     #expressed here using the fact that v_th/w_p = root(2) * Debye length
     #Note that the np.sqrt(2)'s that pop up in relation to the thermal velocities
     #are because some expressions require the RMS velocity
     alpha_e = wpe/(np.sqrt(2)*k*vTe)
     alpha_i = (np.outer(wpi/np.sqrt(2)/vTe, 1/k)).astype(np.cdouble)
     
 
    
     #***********************************
     #Calculate the normalized phase velocities and succeptabilities
     #***********************************
     #See Schaeffer 5.22 for succeptibilities for a Maxwellian plasma
     
     #First for electrons 
     #The second term corrects for drift between the electrons and ions
     # as described in Section 5.3.3 in Sheffield
     xe = wdoppler/(k*vTe) - eidrift/(np.sqrt(2)*vTe)*np.cos(gamma) #Schaeffer 5.21
     
     chiE = -0.5*np.power(alpha_e,2)*ZPrime(xe)#Eq. 5.22 in Schaeffer Thesis
     
     #Then for each species of ions
     xi=1/np.sqrt(2)*np.outer(1/vTi, wdoppler/k) #Schaeffer 5.21

    
     chiI = np.zeros([ion_fract.size, ws.size], dtype=np.cdouble)
     for m in range(ion_fract.size):
          #Eq. 5.22 in Schaeffer Thesis, note that Z Te/Ti is absorbed into alpha_i**2
          chiI[m,:] = -0.5*np.power(alpha_i[m,:],2)*ZPrime(xi[m,:])
          
     #Add the individual ion succeptibilities together 
     chiItot = np.sum(chiI, axis=0)
     epsilon = 1 + chiE + chiItot #Total dielectic function Sheffield Section 5.1
    
     
     #Calculate the spectral density function or "form factor" 
     #Schaeffer Eq. 5.23 and 5.24
     
    
     #Start with just the electron term
     econtr  = 2*np.sqrt(np.pi)/(k*vTe)*np.power(1 - chiE/epsilon,2)*np.exp(-xe**2)
     
     
     
     #Then add on each of the ion terms
     icontr = np.zeros([ion_fract.size, ws.size], dtype=np.cdouble)
     for m in range(ion_fract.size):
         icontr[m,:] = 2*np.sqrt(np.pi)*ion_Z[m]/(k*vTi[m])*np.power(chiE/epsilon,2)*np.exp(-xi[m,:]**2)                                            
     
        
     icontr = np.sum(icontr, axis=0) 
     
     Skw = econtr + icontr
     
     return ks, Skw



def signal(k, spectral_dist, intensity=None, tlen=None, vol=None,  
           probe_wavelength=None, ne=None, solid_angle=4*np.pi):
    """
    Estimate the actual signal as the number of photons detected at a detector
    using experimental parameters.
    
    intensity (float):
        Probe laser intensity in W/cm^2
        
    tlen (float):
        The shorter of either the probe laser pulse length or the scattering timescale in ns
        
    vol (float):
        Scattering volume in cm^3
    
    """
    #Define some constants
    C=2.99792458e10 #Speed of light, cm/s
    thomson_crossection = 6.65e-25 # cm^2
    
    #Unit Conversions
    tlen *= 1e-9 #Convert ns -> s
    probe_wavelength *= 1e-7 #convert nm to cm
    
    #speed of light * plancks constant in J cm^2 * probe wavelength
    hv = 1.98e-23/probe_wavelength

    #Determine the frequency gradient 
    dw = C*np.abs(np.gradient(k)) 
    #Numerically integrate the spectral distrbution function
    int_skw = np.real(np.sum(spectral_dist*dw))
    
    nphotons = tlen/(2*np.pi*hv)*intensity*ne*vol*thomson_crossection*solid_angle*int_skw
    
    return int(nphotons)
     
 
    




"""
print('int S(w): {}'.format(np.sum(FF)))
     
     #Scattered power per solid angle eg. 1.7.13? 
     #Last term (angle correction) is for polarized radiation, Eq. 1.7.14
     #Replace with Eq. 1.7.15 for unpolarized probe radiation
     r = ne*(re**2)*FF*(1 - np.sin(sa)**2*np.cos(dphi)**2)
     
     
     if inst_fcn is not None:
         #Evaluate inst_fcn over a wavelength array chosen to make the output
         #the same length as r (notice that this one must be centered on zero...)
         wspan =  (np.max(wavelength) - np.min(wavelength))/2
         eval_w = np.linspace(-wspan, wspan, num=r.size)
         
         inst_fcn_arr = inst_fcn(eval_w)
         
         #Convolve
         #Note: linear not circular convolution here is important
         r = np.convolve(r, inst_fcn_arr, mode='same')


     #TODO: Roughly estimate the width of the ion feature somehow and block
     #it in a way that doesn't require user input?
     if blockw is not None:
         r = np.where(np.abs(wavelength - laser_wavelength)< blockw, 0, r)

"""



if __name__ == "__main__":
     

  
     ne = 1e15
     
     wavelength = np.arange(520, 545, 0.001)
     
     
     k, spectral_dist = spectral_density(wavelength, Te=10e-3, Ti=1e-3, 
               ion_fract=np.array([1.0]), ion_Z=np.array([1]), 
               ion_Mu=np.array([4]), ne=ne, probe_wavelength=532,
               sa = 63, vpar=0.0, vperp=0.0, eidrift=0.0, gamma=0)
     
     
     
     fig, ax = plt.subplots()

     ax.plot(wavelength, spectral_dist)
     ax.set_xlabel("$\lambda$ (nm)")
     ax.set_ylabel("S(k,w) (s)")
 


     #signal = signal(k, spectral_dist, intensity =1e11, vol=5e-4, solid_angle=0.01, tlen=5, probe_wavelength=532, ne=ne)
     #print("N Photons ~ {:e}".format(signal))
