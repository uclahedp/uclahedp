# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:53:57 2019

@author: Peter
"""

import os
import numpy as np
import matplotlib.pyplot as plt
#from plasmapy.mathematics import plasma_dispersion_func_deriv as ZprimeFcn




def formFactor(Te, fract, Ti, Z, Mu, ne,
                   laser_wavelength, sa, dphi,
                   vpar=0.0, vperp=0.0, eidrift=0.0, gamma=0,
                   wvlgth_range=(400,600,1),
                   blockw = None, inst_fcn=None):
     """
     This function was adapted from a program provided in the book
    "Plasma Scattering of Electromagnetic Radiation" by Froula et al.
     
     Te (float): Electron temperature in KeV
     
     fract (float ndarray): Relative fractions of each ion species
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
     
     wvlgth_range (three float tuple): (start, stop, increment) wavelength
     range in nanometers
     
     
     blockw (None or float): If not none, block a range of width blockw
     around the laser wavelength to zero
     
     inst_fcn (None or function): If not none, represents an insturment
     function (as a python function that takes a wavelength argument and
     returns a normalized insturment function)
     
     """
     
     if np.sum(fract) != 1.0:
          print("WARNING: Sum(IonFractions) != 1.0. Normalizing array." )
          fract = fract/np.sum(fract)
          
    
     #Define Constants
     C=2.99792458e10 #Speed of light, cm/s
     re=2.8179e-13 #Classical Electron radius in cm
     Me=510.9896/C**2 #Electron mass, KeV/C^2
     Mp=Me*1836.1 #Proton mass
     Mi=Mu*Mp #Mass of each ion species
     
     #Define some composite constants that will be used later
     Esq = Me*C**2*re
     const = np.sqrt(4*np.pi*Esq/Me)
     
     #Calculate some plasma values
     Zbar = np.sum(Z*fract)
     ni = fract*ne/Zbar
     wpe = const*np.sqrt(ne) #Electron plasma frequency in Hz
     wpi = const*Z*np.sqrt(ni*Me/Mi)
     vTe = np.sqrt(Te/Me)
     vTi = np.sqrt(Ti/Mi)
     
     
     #Convert all the angles to radians
     sa = sa*2*np.pi/360.0
     dphi = dphi*2*np.pi/360.0
     gamma = gamma*2*np.pi/360.0
     
     #Define the range of scattered frequencies over which we are calculating
     #the spectrum
     wavelength = np.arange(wvlgth_range[0],
                                wvlgth_range[1], 
                                wvlgth_range[2])
     ws = 2*np.pi*C*1e7/wavelength
     
     wl = 2*np.pi*C*1e7/laser_wavelength
     
     #print("{:.2E}, {:.2E}".format(ws[0],ws[-1]))
     

     #Calculate frequencies and wave numbers
     w = ws - wl
     #TODO Explain!
     ks = np.sqrt(ws**2 - wpe**2)/C
     kl = np.sqrt(wl**2 - wpe**2)/C
     #TODO Explain!
     k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(sa))
     #TODO Explain!
     kdotv = (kl - ks*np.cos(sa))*vpar - ks*np.sin(sa)*vperp
     wdoppler = w - kdotv
     
     #TODO WHat are these? Wavenumber normalized to electron and ion inertial
     #lengths?
     klde = (vTe/wpe)*k
     kldi = (np.outer(vTi/wpi, k)).astype(np.cdouble) #TODO is np.outer doing the right thing here?
     
     
     
     #TODO explain this
     #Calculating normalized phase velocities for electrons
     xd = eidrift/(np.sqrt(2)*vTe)*np.cos(gamma)
     xie = wdoppler/(k*np.sqrt(2)*vTe) - xd
     
     #TODO resolve issue with poles in the dispersion function...
     #May need to implement the Cauchy method used in the original after all
     Zpe = ZPrime(xie)
     chiE = -0.5/(klde**2)*Zpe
     

     #Now do the same thing for the ions
     xii=1/np.sqrt(2)*np.outer(1/vTi, wdoppler/k)

     chiI = np.zeros([fract.size, ws.size], dtype=np.cdouble)
     for m in range(fract.size):
          Zpi = ZPrime(xii[m,:])
          chiI[m,:] = -0.5*Zpi/np.power(kldi[m,:],2)
     chiItot = np.sum(chiI, axis=0)
     
     
     #Formfactor
     econtr = (np.sqrt(2*np.pi)/(vTe*k)*
               np.exp(-np.power(xie,2))*
               np.power(np.abs((1 + chiItot)/(1 + chiItot + chiE)),2))
     
     icontr = (2*Ti/Te*np.power(klde,2)/wdoppler*
               np.power(np.abs(chiE),2)/
               np.power(np.abs(1 + chiItot + chiE),2)*np.imag(chiItot))
     

     FF = econtr + icontr
     r = ne*(1 - np.sin(sa)**2*np.cos(dphi)**2)*FF*re**2
     
     
     #TODO figure out what the units are for the form factor??
     
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

     
     return  r



def genZPrimeTable(filepath):
     """
     This function was adapted from a program provided in the book
     "Plasma Scattering of Electromagnetic Radiation" by Froula et al.

     """

     #These values should generally be sufficient
     ximin = -10
     ximax = 10
     dxi = 0.01
     xi = np.arange(ximin, ximax, dxi)
     
     
     N= 4 #Numerical precision
     
     
     #Output array contains three arrays
     #0 -> xi
     #1 -> Real part of Zprime
     #2 -> Imaginary part of Zprime 
     output = np.zeros([3, xi.size])
     output[0,:] = xi
     
     IPV = np.zeros([xi.size])
     RP = np.zeros([xi.size])

     
     for i in range(xi.size):
          #phi is min distance from a singularity
          phi = dxi*np.abs(xi[i]) + 1e-6
          dz = phi/N
          
          #Make arrays of values to calculate on either side of xi
          zm = np.arange(xi[i] - phi, ximin-1, -dz)
          zp = np.arange(xi[i] + phi, ximax+1, dz)
          
          #Perform the Riemann sum integral over each of these ranges
          sum_zp = dz*np.sum(zp*np.exp(-np.power(zp,2))/(zp-xi[i]))   
          sum_zm = dz*np.sum(zm*np.exp(-np.power(zm,2))/(zm-xi[i]))
          
          IPV[i] = sum_zp + sum_zm
          RP[i] = 2*phi*(1-2*np.power(xi[i],2))
          
     dW = 1/np.sqrt(np.pi)*(IPV + np.exp(-np.power(xi,2))*(RP-1j*np.pi*xi))
     
     output[1,:] = -2*np.real(dW)
     output[2,:] = 2*np.imag(dW)
     
     np.savetxt(filepath, output)
     
   


def ZPrime(xi):
    
     """
     This function was adapted from a program provided in the book
     "Plasma Scattering of Electromagnetic Radiation" by Froula et al.
    
     If the file zprime.txt does not exist in the current directory, generate it

    
     """
     
     #This call gets the directory of THIS script, not the cwd
     #so that the zprime.txt file lives in UCLAHEDP's folder, not wherever
     #the user is running their own code.
     scriptdir = os.path.dirname( os.path.realpath(__file__) )

     table_path = os.path.join(scriptdir, "zprime.txt")
     if not os.path.exists(table_path):
         print("Zprime.txt file not found: generating")
         genZPrimeTable(table_path)
          
     table = np.loadtxt(table_path)
     
     xitable = table[0,:]
     
     rZptable = table[1,:]
     iZptable = table[2,:]
     
     #plt.plot(xitable, rZptable)
     #plt.show()
     
     
     ai = np.argmin(np.abs(xi - np.min(xitable)))
     bi = np.argmin(np.abs(xi - np.max(xitable)))

     
     rZp = np.zeros([xi.size])
     iZp = np.zeros([xi.size])
     
     
     #xi range below table range
     rZp[0:ai] = np.power(xi[0:ai],-2)
     iZp[0:ai] = 0.0
     
     #xi range above table range
     rZp[bi:-1] = np.power(xi[bi:-1],-2)
     iZp[bi:-1] = 0.0
     
     
     #xi range within table range
     rZp[ai:bi] = np.interp(xi[ai:bi], xitable, rZptable)
     iZp[ai:bi] = np.interp(xi[ai:bi], xitable, iZptable)
     
     return rZp + 1j*iZp
     



     

if __name__ == "__main__":
     

     Te = 5e-3
     fract= np.array([1.0])
     Ti = 1e-3
     Z = np.array([1])
     Mu = np.array([4])
     
     laser_wavelength = 532 #3.5e15 rad/s 532 nm
     
     wvlgth_range=(500,600,.01)
     
     wavelength, r = formFactor(Te, fract , Ti, Z, Mu, 1e13,
                   laser_wavelength, 63, 90, wvlgth_range=wvlgth_range,
                   vpar=0, vperp=0, eidrift=0, gamma=0)
     
     r = np.where(np.abs(wavelength - laser_wavelength)< 0.2, 0, r)

     plt.plot(wavelength, r)
     plt.xlim(525,540)
     #plt.ylim(0,1e-24)
