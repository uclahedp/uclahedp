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
from scipy.special import iv as BesselIV


def spectrum(wavelength, mode='collective', Te=None, Ti=None, 
                   ion_fract=np.array([1.0]), ion_Z=None, 
                   ion_Mu=None, ne=None, probe_wavelength=None, B0 = np.array([0,0,0]),
                   probe_n = np.array([1,0,0]), scatter_n = np.array([0,1,0]), pol_n=None,
                   ion_v = np.array([0,0,0]), electron_v = np.array([0,0,0]),
                   scattering_length = None, scattered_power = False, inst_fcn=None,
                   block_width=None, verbose=False):
     
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
         Plasma mean  (0th order) electron number density in cm^-3. 
     
     probe_wavelength (float): 
         Incident probe laser wavelength in nm
         
     B0 (float ndarray):
         Vector magnetic field in Gauss. Default value is (0,0,0) which is interpreted
         as unmagnetized scattering.
         
     probe_n (float ndarray):
         Normalized unit vector representing the orientation of the probe laser beam
         
    
     scatter_n (float ndarray):
         Normalized unit vector representing the path of scattered light to the detector
         
     pol_n (float ndarray):
         Normalized unit vector representing the polarization of the probe laser beam. Setting
         this variable to None (the default) represents an unpolarized beam.
         
     
     ion_v (float ndarray):
         Velocity vector  in cm/s for the ion population
         #TODO: support multiple ion population velocities for different species
         
     electron_v (float ndarray):
         Velocity vector in cm/s for the electron population
         
         
     scattering_length (float):
         Length of the scattering region in cm
         
         
     scattered_power (boolean):
         If set to True, the function returns the scattered power spectrum
         normalized to the incident power : P_s/P_i. Note this calculation requires
         the scattering_length variable to be set.
         (Eq. 5.15 in Schaeffer or 5.1.1 in Sheffield)
         
         
     inst_fcn (function):
         An insturment function to apply to the spectral density function to simulate the spectrum
         as it would be observed by a detector.
         
         
     block_width (float):
         Set a range of [-block_width, block_width] to zero in the spectral density function.
         This is useful for blocking out the ion feature when fitting the electron feature.
   
     """
     
     #TODO: throw an error if the lengths of any of the arrays are inconsistent
     #with the number of ion species (i.e Z, mu, V, Ti)
     

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
     
     #k is the wavenumber difference (as illustrated in Fig.1.5 of Sheffield)
     sa = np.arccos(np.dot(probe_n, scatter_n)) #rad
     k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(sa)) #Eq. 1.7.10 in Sheffield, units 1/cm
     k_n = scatter_n - probe_n #Normal vector along k
     
     #Compute Doppler-shifted frequencies for both the ions and electrons
     #TODO: Is specifying both population's velocities in the rest frame
     #the best choice?
     wdoppler_i = w - np.dot( np.outer(k, k_n), ion_v)
     wdoppler_e = w - np.dot( np.outer(k, k_n), electron_v)
     
         
     #Compute the scattering parameter alpha
     #expressed here using the fact that v_th/w_p = root(2) * Debye length
     #Note that the np.sqrt(2)'s that pop up in relation to the thermal velocities
     #are because only some expressions require the RMS velocity.
     alpha = wpe/(np.sqrt(2)*k*vTe)
     
     
     if verbose:
         print("min/max alpha: {:.3f}, {:.3f}".format(np.min(alpha),np.max(alpha)))


     #************************************************************************
     #Non-Collective Scattering, Unmagnetized 
     #************************************************************************
     
     #TODO: include a non-collective magnetized case based on Sheffield Sec. 4.6
     
     if mode == 'non-collective':
         xe = wdoppler_e/(k*np.sqrt(2)*vTe)
         #Schaeffer Eq. 5.26, Sheffield Eq. 4.5.1 (PDF pg. 84)
         Skw = 2*np.pi*np.exp(-xe**2)/(np.sqrt(np.pi)*k*np.sqrt(2)*vTe)
     
     
        
        
     #************************************************************************
     #Collective Scattering, Unmagnetized 
     #************************************************************************
        
     elif bmag == 0:
         #Calculate the normalized phase velocities and succeptabilities
         #Following section Sec. 5.1.4 in Schaeffer's thesis
         #and Sec. 3.4.2 in Sheffield
         xe = wdoppler_e/(k*np.sqrt(2)*vTe)
         xi=1/np.sqrt(2)*np.outer(1/vTi, wdoppler_i/k) #Schaeffer 5.21

         
         if verbose:
             print("min/max x_e: {:.3f}, {:.3f}".format(np.min(xe),np.max(xe)))
             print("min/max x_i: {:.3f}, {:.3f}".format(np.min(xi),np.max(xi)))
         
         chiE = -0.5*np.power(alpha,2)*ZPrime(xe)#Eq. 5.22 in Schaeffer Thesis
         
         
         #Treatment of multiple ion species as described in Sheffield Sec. 5.1
         chiI= np.zeros([ion_fract.size, w.size], dtype=np.cdouble)
         for m in range(ion_fract.size):
              #Eq. 5.22 in Schaeffer Thesis
              chiI[m,:] = -0.5*np.power(alpha,2)*ion_Z[m]*Te/Ti*ZPrime(xi[m,:])
         
         epsilon = 1 + chiE + np.sum(chiI, axis=0)
        
        
         #Calculation of the spectral dispersion function
         #Schaeffer 
         #or Sheffield Sec. 5.1
         econtr = 2*np.sqrt(np.pi)/k/vTe*np.power(np.abs(1 - chiE/epsilon),2)*np.exp(-xe**2)
        
         icontr = np.zeros([ion_fract.size, w.size], dtype=np.cdouble) 
         for m in range(ion_fract.size):  
            icontr[m,:] = 2*np.sqrt(np.pi)*ion_Z[m]/k/vTi*np.power(np.abs(chiE/epsilon),2)*np.exp(-xi[m,:]**2)
         
         #Re-cast as real as a formality: imaginary part is zero
         Skw = np.real(econtr + np.sum(icontr, axis=0))
         

     #************************************************************************
     #Collective Scattering, Magnetized Electrons and Ions
     #************************************************************************
        
     elif bmag > 0: 
        print("Magnetized Thomson Scattering is only roughly implemented: use at own risk")
        
        #TODO: finish this section. This calculation is numerically more complicated than
        #the unmagnetized case.
         
        #Compute vperp/vpara and gyroradii
        k_para = k*np.dot(k_n, b_n)
        k_perp = np.sqrt(k**2 - k_para**2)
        rho_e = vTe/wce
        rho_i = np.mean(vTi/wci)
        
        
        
        larr = np.arange(-10, 10) #Array of resonance numbers to include in calculation

        earg = k_perp**2*rho_e**2
        iarg = k_perp**2*rho_i**2
        
        
        #TODO: When iarg is too big (corresponding to unmagnetized ions)
        #switch models to avoid numerical problems with large args in Bessell fcns.
        
        if verbose:
             print("min/max earg: {:.3f}, {:.3f}".format(np.min(earg),np.max(earg)))
             print("min/max iarg: {:.3f}, {:.3f}".format(np.min(iarg),np.max(iarg)))
        
    
        #He, Sheffield 10.3.5
        He = np.complex128(alpha**2)
        for l in larr:
            t1 = alpha**2*np.exp(-earg)*BesselIV(l,earg)
            t2 = wdoppler_e / (wdoppler_e - l*wce)
            xel = (wdoppler_e - l*wce)/(k_para*np.sqrt(2)*vTe)
            t3 = -ZFcn(xel)*xel
            He += -t1*t2*t3
            
        
        #Hi, Sheffield 10.3.7
        Hi = np.complex128(alpha**2*mean_Z*Te/Ti)
        for l in larr:
            t1 = alpha**2*mean_Z*Te/Ti*np.exp(-iarg)*BesselIV(l,iarg)
            t2 = wdoppler_i / (wdoppler_i - l*wci)
            xil = (wdoppler_i - l*wci)/(k_para*np.sqrt(2)*vTi)
            t3 = -ZFcn(xil)*xil
            Hi += -t1*t2*t3
            
        epsilon = 1 + He + Hi
            
        
        #Compute Skw, Sheffield 10.3.9
        sum_e = np.complex128(0)    
        sum_i = np.complex128(0)
        for l in larr:
            xel = (wdoppler_e - l*wce)/(k_para*np.sqrt(2)*vTe)
            t1 = np.exp(-earg)*BesselIV(l,earg)
            t2 = np.exp(-xel**2)/k_para/(np.sqrt(2)*vTe)
            sum_e += t1*t2
            
            xil = (wdoppler_i - l*wci)/(k_para*np.sqrt(2)*vTi)
            t1 = np.exp(-iarg)*BesselIV(l,iarg)
            t2 = np.exp(-xil**2)/k_para/(np.sqrt(2)*vTi)
            sum_i += t1*t2
            
            
        econtr = 2*np.sqrt(np.pi)*np.power(np.abs(1 - He/epsilon),2)*sum_e
        
        icontr = 2*np.sqrt(np.pi)*mean_Z*np.power(np.abs(He/epsilon),2)*sum_i
        
        Skw = np.real(econtr + icontr)
        
        
        
        
        
     #************************************************************************
     #Other final options based on keywords
     #************************************************************************
     if inst_fcn is not None:
         #Evaluate inst_fcn over a wavelength array chosen to make the output
         #the same length as r (notice that this one must be centered on zero...)
         wspan =  (np.max(wavelength) - np.min(wavelength))/2
         eval_w = np.linspace(-wspan, wspan, num=Skw.size)
         
         inst_fcn_arr = inst_fcn(eval_w)
         
         #Convolve
         #Note: linear not circular convolution here is important
         Skw = np.convolve(Skw, inst_fcn_arr, mode='same')
      
        
     if block_width is not None:
         Skw = np.where(np.abs(wavelength - probe_wavelength)< block_width, 0, Skw) 
        
        
        
     if scattered_power:
        
        if pol_n is not None:
            
            #Create a vector perpendicular to the scattering plane
            temp = np.cross(scatter_n, probe_n)
            temp = temp/np.linalg.norm(temp) #normalize

            #Compute dphi as defined in Fig. 1.6 of Sheffield
            dphi = np.pi/2  - np.arccos( np.dot(pol_n, temp) )  
            
            #Sheffield Eq. 1.7.14
            pterm = 1 - np.sin(sa)**2*np.cos(dphi)**2
        
        else:
            #Sheffield Eq. 1.7.15
            pterm = 1 - 0.5*np.sin(sa)**2
            
        #Compute the "geometric factor"
        gf = 1 + 2*wdoppler_e/wl
            
        #Compute the scattered power
        #(Eq. 5.15 in Schaeffer or 5.1.1 in Sheffield)
        PS = re**2*gf*scattering_length*ne/(2*np.pi)*pterm**2*Skw
        return ks, PS
    
    
     else:
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
     
 
    
def _gaussian_inst_fcn(wavelength):
    #Create an instrument function that represents the spectral response
    #of the spectrometer recording the Thomson spectrum
    res= 1
    fcn = np.exp(-0.5*np.power(wavelength/res,2))/(res*np.sqrt(2*np.pi))
    return fcn    






if __name__ == "__main__":

     ne = 1e13
     
     
     w = 532
     pmw = 100
     #wavelength = np.arange(w-pmw, w+pmw, 0.01)
     wavelength = np.arange(520,545, 0.01)
     
     sa = 90
     probe_n = np.array([1,0,0])
     scatter_n = np.array([np.cos(np.deg2rad(sa)), np.sin(np.deg2rad(sa)), 0])
     
     
     pol_n = np.array([0,0,1])
     
     electron_v = np.array([0,0,0])
     ion_v = np.array([0,0,0])
     
     ba = 70
     bmag = 0#7e6
     B0 = np.array([0,bmag*np.cos(np.deg2rad(ba)),bmag*np.sin(np.deg2rad(ba))])
     
     
     k, spectral_dist = spectrum(wavelength, mode='collective',
               Te=10e-3, Ti=1e-3, electron_v=electron_v, ion_v = ion_v,
               probe_n = probe_n, scatter_n=scatter_n, B0=B0,
               pol_n = None,
               ion_fract=np.array([1.0]), ion_Z=np.array([1]), 
               ion_Mu=np.array([1]), ne=ne, probe_wavelength=w,
               scattered_power=False, scattering_length=0.5, 
               block_width  = None, inst_fcn = None, verbose=True)
     
     
     
     fig, ax = plt.subplots()
     ax.set_xlabel("$\lambda$ (nm)")
     ax.set_ylabel("S(k,w)")
     ax.axvline(x=w, color='red')
     
     
     ax.set_ylim(0, 1.2*np.max(spectral_dist))
     #ax.set_ylim(0, .5e-13)
     ax.set_xlim(np.min(wavelength), np.max(wavelength))

     ax.plot(wavelength, spectral_dist, lw=1)
     
     
 


     signal = signal_estimate(k, spectral_dist, intensity =1e11, vol=5e-4, solid_angle=0.01, tlen=5, probe_wavelength=532, ne=ne)
     print("N Photons ~ {:.1e}".format(signal))
