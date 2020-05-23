# -*- coding: utf-8 -*-
"""
thomson_scattering.py 

@author: Peter Heuer

This module contains code for the analysis of Thomson scattering spectra
from plasmas with multiple ion species.

This function was adapted from a program provided in the book
     "Plasma Scattering of Electromagnetic Radiation" by Froula et al.

"""

import numpy as np
import matplotlib.pyplot as plt
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func_deriv as ZPrime
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func as ZFcn

from scipy.special import iv as BesselI


def spectral_density(wavelength, mode='collective', Te=None, Ti=None, 
                   ion_fract=np.array([1.0]), ion_Z=None, 
                   ion_Mu=None, ne=None, probe_wavelength=None, B0 = np.array([0,0,0]),
                   probe_n = np.array([1,0,0]), scatter_n = np.array([0,1,0]), 
                   ion_v = np.array([0,0,0]), electron_v = np.array([0,0,0]),
                   edrift = np.array([0,0,0]) ):
     
     """
     model (str):
         Select the model to be used for S(k,w). Choices are:
             "non-collective"
             "collective"
     
     
     Te (float): 
         Electron temperature in KeV
     
     Ti (float): 
         Ion temperature in KeV (all species must have the same temp?)
     
     ion_fract (float ndarray): 
         Relative fractions of each ion species present. Must sum to 1.0.
     
     ion_Z (float ndarray): 
         Ion charge states, in units of the proton charge
     
     ion_Mu (float ndarray): 
         Ion atomic masses in units of the proton mass
     
     ne (float): 
         Plasma electron number density in cm^-3. 
     
     probe_wavelength (float): 
         Incident probe laser wavelength in nm
         
     
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
     

     
     """
     
     #TODO: throw an error if the lengths of any of the arrays are inconsistent
     #with the number of ion species (i.e Z, mu, V, Ti)
     
     #TODO: Support multiple ion populations with different vpar and vperp
     

     
     if np.sum(ion_fract) != 1.0:
          print("WARNING: Sum(IonFractions) != 1.0. Normalizing array." )
          ion_fract = ion_fract/np.sum(ion_fract)
          
          
     #Make sure all the unit vectors are normalized
     probe_n = probe_n/np.linalg.norm(probe_n)
     scatter_n = scatter_n/np.linalg.norm(scatter_n)
     
     bmag = np.linalg.norm(B0)
     if bmag > 0:
         b_n = B0/bmag
       
    
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
     mean_Mu = np.sum(ion_Mu*ion_fract)
     ni = ion_fract*ne/mean_Z #cm^-3
     wpe = const*np.sqrt(ne) #Electron plasma frequency in rad/s
     wpi = const*ion_Z*np.sqrt(ni*Me/Mi) # rad/s
     wce = 1.7e7*bmag #rad/s
     wci = 9.5e3*mean_Z/mean_Mu*bmag #rad/s
     vTe = np.sqrt(Te/Me) #cm/s
     vTi = np.sqrt(Ti/Mi) #cm/s
     
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
     sa = np.arccos(np.dot(probe_n, scatter_n))
     k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(sa)) #Eq. 1.7.10 in Sheffield, units 1/cm
     k_n = scatter_n - probe_n #Normal vector along k
     
     #Compute Doppler-shifted frequencies for both the ions and electrons
     #TODO: Is specifying both population's velocities in the rest frame
     #the best choice?
     #TODO: In future support different fluid velocities for each ion species?
     wdoppler_i = w - np.dot( np.outer(k, k_n), ion_v)
     wdoppler_e = w - np.dot( np.outer(k, k_n), electron_v)
         
     #Compute the scattering parameter alpha
     #expressed here using the fact that v_th/w_p = root(2) * Debye length
     #Note that the np.sqrt(2)'s that pop up in relation to the thermal velocities
     #are because only some expressions require the RMS velocity.
     alpha_e = wpe/(np.sqrt(2)*k*vTe)
     alpha_i = (np.outer(wpi/np.sqrt(2)/vTe, 1/k)).astype(np.cdouble)
     

     #************************************************************************
     #Non-Collective Scattering, Unmagnetized 
     #************************************************************************
     
     #TODO: include a non-collective, magnetized case?
     
     if mode == 'non-collective':
         #Schaeffer Eq. 5.26
         Skw = 2*np.pi*np.exp(-np.power(wdoppler_e/k/vTe,2))/(np.sqrt(np.pi)*k*vTe)
         return ks, Skw
     
     
     #************************************************************************
     #Collective Scattering, Unmagnetized 
     #************************************************************************
        
     if bmag == 0:
         #Calculate the normalized phase velocities and succeptabilities
         #See Schaeffer 5.22 for succeptibilities for a Maxwellian plasma
     
         #First for electrons 
         #The second term corrects for drift between the electrons and ions
         # as described in Section 5.3.3 in Sheffield and #Schaeffer Eq. 5.21
         xe = wdoppler_e/(k*vTe) 
         chiE = -0.5*np.power(alpha_e,2)*ZPrime(xe)#Eq. 5.22 in Schaeffer Thesis
         
         
         #Then for each species of ions
         xi=1/np.sqrt(2)*np.outer(1/vTi, wdoppler_i/k) #Schaeffer 5.21
    
         chiI = np.zeros([ion_fract.size, ws.size], dtype=np.cdouble)
         for m in range(ion_fract.size):
              #Eq. 5.22 in Schaeffer Thesis, note that Z Te/Ti is absorbed into alpha_i**2
              chiI[m,:] = -0.5*np.power(alpha_i[m,:],2)*ZPrime(xi[m,:])
              
         #Add the individual ion succeptibilities together 
         chiItot = np.sum(chiI, axis=0)
         epsilon = 1 + chiE + chiItot #Total dielectic function Sheffield Section 5.1
         
    
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
     
     #************************************************************************
     #Collective Scattering, Magnetized
     #************************************************************************
        
     elif bmag > 0: 
        #Compute vperp/vpara and gyroradii
        k_para = k*np.dot(k_n, b_n)
        k_perp = np.sqrt(k**2 - k_para**2)
        rho_e = vTe/wce/np.sqrt(2)
        rho_i = vTi/wci/np.sqrt(2)
        
        larr = np.arange(-3, 3) #Array of resonance numbers to include in calculation
    
        earg = k_perp**2*rho_e**2
        iarg = k_perp**2*rho_i**2
        
 
    
        #First compute He
        const = alpha_e**2
        He = const
        for l in larr:     
            term = const*np.exp(-earg)*BesselI(l,earg)*wdoppler_e/(wdoppler_e - l*wce)
       
            if any(k_para != 0):
                xel = (wdoppler_e - l*wce)/(k_para*np.sqrt(2)*vTe)
                #Term in brackets in Sheffield Eq. 10.3.6 is basically the plasma disp. fcn.
                term *= -xel*ZFcn(xel)
            He += -term
            
        #Next compute Hi
        const = mean_Z*Te/Ti*np.mean(alpha_i,axis=0)**2
        Hi = const
        
        for l in larr:
            term = const*np.exp(-iarg)*BesselI(l,iarg)*wdoppler_i/(wdoppler_i - l*wci)
            if any(k_para != 0):
                xil = (wdoppler_i - l*wci)/(k_para*np.sqrt(2)*vTi)
                term *= -xil*ZFcn(xil)
            Hi += -term
                 
        #Now the total epsilon
        epsilon = 1 + He + Hi
        
        #Now compute the spectral dispersion function
        esum, isum = 0,0
        for l in larr:
            earg2 = (wdoppler_e - l*wce)/(k_para*vTe)
            iarg2 = (wdoppler_i - l*wci)/(k_para*vTi)
            
            esum += np.exp(earg)*BesselI(l,earg)*np.exp(-earg2**2)/(k_para*vTe)
            isum += np.exp(iarg)*BesselI(l,iarg)*np.exp(-iarg2**2)/(k_para*vTi)
           
        
        Skw = ( 2*np.sqrt(np.pi)*np.power(1 - He/epsilon,2)*esum + 
                2*np.sqrt(np.pi)*mean_Z*np.power(He/epsilon,2)*isum )
               
        print(Skw)
        
        return ks, Skw
                
        
        
         
    



def signal_estimate(k, spectral_dist, intensity=None, tlen=None, vol=None,  
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
blockw (None or float): If not none, block a range of width blockw
     around the laser wavelength to zero
     
     inst_fcn (None or function): If not none, represents an insturment
     function (as a python function that takes a wavelength argument and
     returns a normalized insturment function)

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
     
     
     w = 532
     wavelength = np.arange(w-10, w+10, 0.001)
     
     sa = 63
     probe_n = np.array([1,0,0])
     scatter_n = np.array([np.cos(np.deg2rad(sa)), np.sin(np.deg2rad(sa)), 0])
     
     electron_v = np.array([0,0,0])
     ion_v = np.array([0,0,0])
     
     B0 = np.array([0,0,0])
     
     
     k, spectral_dist = spectral_density(wavelength, mode='collective',
               Te=5e-3, Ti=1e-3, electron_v=electron_v, ion_v = ion_v,
               probe_n = probe_n, scatter_n=scatter_n, B0=B0,
               ion_fract=np.array([1.0]), ion_Z=np.array([1]), 
               ion_Mu=np.array([4]), ne=ne, probe_wavelength=w)
     
     
     
     fig, ax = plt.subplots()

     ax.plot(wavelength, np.real(spectral_dist))
     ax.set_xlabel("$\lambda$ (nm)")
     ax.set_ylabel("S(k,w) (s)")
 


     signal = signal_estimate(k, spectral_dist, intensity =1e11, vol=5e-4, solid_angle=0.01, tlen=5, probe_wavelength=532, ne=ne)
     print("N Photons ~ {:.1e}".format(signal))
