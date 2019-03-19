# -*- coding: utf-8 -*-
"""
Generates theoretical dispersion relation curves for parallel beam instabilities
@author: Peter Heuer
"""
import os
import numpy as np
import scipy.optimize as optimize
import scipy.interpolate as interpolate
from scipy.signal import savgol_filter, medfilt
from scipy.signal import find_peaks

import matplotlib.pyplot as plt

def vb_to_vbLab(nb, vb, qb, qc):
     vc = -(qb/qc)*nb*vb/(1-qb*nb)
     vlab = np.abs(vc - vb)
     return vlab
 
    
def vbLab_to_vb(nb, vblab, qb, qc):
    vb = np.abs( vblab / (1 + qb*nb/(qc*(1 - nb)) ) )
    return vb


def dispRel(z,k, vb=None, nb=None, qb=None, qc=None, mb=None):
    if vb is None:
        vb = 5 
    if nb is None:
        nb = 0.1
    if qb is None:
        qb = 4
    if qc is None:
        qc = 1
    if mb is None:
        mb = 3.0
        
    #Construct a complex value
    w = z[0] + 1j*abs(z[1])

    #init output array
    F = np.zeros((2))

    ne = 1.00
    nc = ne - qb*nb
    
    ve = 0.0
    vc = -(qb/qc)*nb*vb/(1-qb*nb)
    
    wcc = 1
    wcb = (qb/mb)*wcc
    #wce -> inf
    
    wpc = 2900.0

    #Calculate solution
    sol = ( (pow(w,2)/pow(wpc,2) - pow(k,2))*(w - vc*k + wcc)*(w - vb*k + wcb) - 
        nc*(w - vc*k)*(w - vb*k + wcb) + 
        ne*( (w - ve*k)/wcc )*(w - vc*k + wcc)*(w - vb*k + wcb) -
        nb*(pow(qb,2)/mb)*(w-vb*k)*(w - vc*k + wcc) )
    
    #Repackage as two reals
    F[0] = np.real(sol)
    F[1] = np.imag(sol)
     
    return F

def bestGuess(k):
    return np.array([k**2, 1])



def solveDispRel(nb=None, vb=None, qb=None, qc=None, mb=None, krange = None, guessFcn = None):
    
    #Create kspace and output w vectors
    if krange is None:
        krange = (-10, 10)
        

    if guessFcn is None:
        guessFcn = bestGuess
        
    kspace = np.arange(krange[0], krange[1],  .01)
    w = np.empty([len(kspace)], dtype=np.complex128)

    #Loop through, solve for each k
    for i,k in enumerate(kspace):
        #Guess: k**2. This is important to get the right branch
        guess = guessFcn(k)
        wr, wi = optimize.fsolve(lambda z : dispRel(z, k, vb=vb, nb=nb, qb=qb, qc=qc, mb=mb), guess, factor=0.1)
        #Repackage as a complex value
        #The absolute value on wi seems to be necessary to keep on the right
        #branch...
        w[i] = wr + 1j*abs(wi)
    
    #Interpolate roots
    k = np.arange(krange[0], krange[1], .01)
    wr = interpolate.interp1d(kspace, np.real(w) ) (k)
    wi = interpolate.interp1d(kspace, np.imag(w) ) (k)
    
    return k, wr, wi


def classify_peaks(k, wi):
    
    dk = np.mean(np.gradient(k))
    kwidth = .01
    width = int(kwidth/dk)
    
    wi = medfilt(wi, 9)
    
    peaks = find_peaks(wi, height = 0.001)
    peaks = peaks[0]

    output  = {'F': {'k':None, 'gr':0},
               'A': {'k':None, 'gr':0},
               'C': {'k':None, 'gr':0},
               'E': {'k':None, 'gr':0}}

    pos = [p for p in peaks if k[p]> 0]
    neg = [p for p in peaks if k[p]< 0]

    if pos == []:
        pass
    elif len(pos) == 1:
        output['F'] = {'k': k[pos[0]], 'gr': wi[pos[0]]}
    elif len(pos) == 2:
        output['F'] = {'k': k[pos[0]], 'gr': wi[pos[0]]}
        output['E'] = {'k': k[pos[1]], 'gr': wi[pos[1]]}
    else:
        raise ValueError("Invalid number of positive k peaks!")
     
    if neg == []:
        pass
    elif len(neg) == 1:
        output['A'] = {'k': k[neg[0]], 'gr': wi[neg[0]]}
    elif len(neg) == 2:
        output['A'] = {'k': k[neg[1]], 'gr': wi[neg[1]]}
        output['C'] = {'k': k[neg[0]], 'gr': wi[neg[0]]}
    else:
        raise ValueError ("Invalid number of negative k peaks!")

    return output




if __name__ == '__main__':
    
    vb = 5
    k, wr, wi =  solveDispRel(nb=0.05, vb=vb, qb=4, krange=(-5,5))
   
    fig, ax = plt.subplots( figsize = [4,4])
    
    vbline = k*vb
    
    #Filtering bc branch in solution makes an annoying bump
    wr = savgol_filter(wr, 501, 3)
    
    wifactor = 25
    ax.plot(k, wr, '-', k, wi*wifactor, '--')
    ax.plot(k, vbline, '--')
    ax.axvline(0, color='black', linewidth=1)
    plt.show()