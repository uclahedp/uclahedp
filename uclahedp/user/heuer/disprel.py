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

from uclahedp.tools.util import timeRemaining as timer

import h5py
import matplotlib.pyplot as plt

def vb_to_vbLab(nb, vb, qb, qc):
     vc = -(qb/qc)*nb*vb/(1-qb*nb)
     vlab = np.abs(vc - vb)
     return vlab
 
    
def vbLab_to_vb(nb, vblab, qb, qc):
    vb = np.abs( vblab / (1 + qb*nb/(qc*(1 - nb)) ) )
    return vb


def calc_vc(nb, vb, qb, qc):
    return (qb/qc)*nb*vb/(1-qb*nb)

def calc_nc(nb, qb):
     return 1.00 - qb*nb


def dispRel(z,k, vb=None, nb=None, qb=None, qc=None, mb=None, mc=None):
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
    if mc is None:
        mc = 1.0
        
    #Construct a complex value
    w = z[0] + 1j*abs(z[1])

    #init output array
    F = np.zeros((2))

    ne = 1.00
    nc = calc_nc(nb, qb)
    
    ve = 0.0
    vc = calc_vc(nb, vb, qb, qc)
    
    wcc = 1
    wcb = (qb/mb)*wcc
    #wce -> inf
    
    wpc = 2900.0

    #Calculate solution
    sol = ( (pow(w,2)/pow(wpc,2) - pow(k,2))*(w - vc*k + wcc)*(w - vb*k + wcb) - 
        nc*(pow(qc,2)/mc)*(w - vc*k)*(w - vb*k + wcb) + 
        ne*( (w - ve*k)/wcc )*(w - vc*k + wcc)*(w - vb*k + wcb) -
        nb*(pow(qb,2)/mb)*(w-vb*k)*(w - vc*k + wcc) )
    
    #Repackage as two reals
    F[0] = np.real(sol)
    F[1] = np.imag(sol)
     
    return F

def bestGuess(k):
    return np.array([k**2, 1])

def rightPropGuess(k):
    if k < 0:
        return np.array([k*0.1, -1])
    else: 
        return np.array([k**2, .1])

def leftPropGuess(k):
    return rightPropGuess(-k)



def solveDispRel(nb=None, vb=None, qb=None, qc=None, mb=None, mc=None, krange = None, guessFcn = None):
    
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
        wr, wi = optimize.fsolve(lambda z : dispRel(z, k, vb=vb, nb=nb, qb=qb, qc=qc, mb=mb, mc=mc), guess, factor=0.1)
        #Repackage as a complex value
        #The absolute value on wi seems to be necessary to keep on the right
        #branch...
        w[i] = wr + 1j*abs(wi)
    
    #Interpolate roots
    k = np.arange(krange[0], krange[1], .01)
    wr = interpolate.interp1d(kspace, np.real(w) ) (k)
    wi = interpolate.interp1d(kspace, np.imag(w) ) (k)
    
    return k, wr, wi


def classify_peaks(k, wr, wi):
    
    dk = np.mean(np.gradient(k))
    kwidth = .01
    width = int(kwidth/dk)
    
    wi = medfilt(wi, 25)
    
    peaks = find_peaks(wi, height = 0.001)
    peaks = peaks[0]

    output  = {'F': {'k':None, 'w':None, 'gr':0},
               'A': {'k':None, 'w':None, 'gr':0},
               'C': {'k':None, 'w':None, 'gr':0},
               'E': {'k':None, 'w':None, 'gr':0}}

    pos = [p for p in peaks if k[p]> 0]
    neg = [p for p in peaks if k[p]< 0]

    if pos == []:
        pass
    elif len(pos) == 1:
        output['F'] = {'k': k[pos[0]], 'w': wr[pos[0]], 'gr': wi[pos[0]]}
    elif len(pos) == 2:
        output['F'] = {'k': k[pos[0]], 'w': wr[pos[0]], 'gr': wi[pos[0]]}
        output['E'] = {'k': k[pos[1]], 'w': wr[pos[1]], 'gr': wi[pos[1]]}
    else:
        print(len(pos))
        raise ValueError("Invalid number of positive k peaks!")
     
    if neg == []:
        pass
    elif len(neg) == 1:
        output['A'] = {'k': k[neg[0]], 'w': wr[neg[0]], 'gr': wi[neg[0]]}
    elif len(neg) == 2:
        output['A'] = {'k': k[neg[1]], 'w': wr[neg[1]], 'gr': wi[neg[1]]}
        output['C'] = {'k': k[neg[0]], 'w': wr[neg[0]], 'gr': wi[neg[0]]}
    else:
        print(len(neg))
        raise ValueError ("Invalid number of negative k peaks!")

    return output



def gen_gr():
    qb = 4.0
    mb = 3.0
    qc = 1
    mc = 1
    
    #NOTE: Avoid letting nb/n0 = 1/qb = 0.25! Jump's branches there...
    
    krange = (-10,10)
    vblab_range = np.arange(1, 20, 0.1)
    nb_range = np.arange(0, 0.24, 0.01)
    
    #vblab_range = np.arange(4, 8, 1)
    #nb_range = np.arange(.1, 0.2, .05 )
    
    nv = len(vblab_range)
    nn = len(nb_range)
     
    shape = (nv, nn)

    output = {}
    output['F'] = {'k':np.zeros(shape), 'w':np.zeros(shape), 'gr':np.zeros(shape)}
    output['A'] = {'k':np.zeros(shape), 'w':np.zeros(shape), 'gr':np.zeros(shape)}
    output['C'] = {'k':np.zeros(shape), 'w':np.zeros(shape), 'gr':np.zeros(shape)}
    output['E'] = {'k':np.zeros(shape), 'w':np.zeros(shape), 'gr':np.zeros(shape)}

    output['vblab'] = vblab_range
    output['nb'] = nb_range
    output['vbprf'] = np.zeros(shape)
    output['vcprf'] = np.zeros(shape)
    
    
    t = timer(nv*nn, reportevery=1)
    
    for i, vblab in enumerate(vblab_range):
        for j, nb in enumerate(nb_range):
            ii = i*nn + j
            t.updateTimeRemaining(ii)
            
            vb  = vbLab_to_vb(nb, vblab, qb, qc)
            vc = calc_vc(nb, vb, qb, qc)
            
            output['vbprf'][i,j] = vb
            output['vcprf'][i,j] = vc
            

            print("i,j=" + str(i) + ','  + str(j) + 
                    ", nb=" + str(np.round(nb, decimals=2)) +
                  ', vblab=' + str(np.round(vblab, decimals=2)) +
                  ', vb=' + str(np.round(vb, decimals=2))  )
            
            k, wr, wi = solveDispRel( nb=nb, vb =vb, qb=qb, mb=mb, qc=qc, mc=mc, guessFcn = bestGuess, krange=krange)
            peaks = classify_peaks(k, wr, wi)
            #plot1(k, wr, wi, peaks=peaks)
            for g in ('F', 'A', 'C', 'E'):
                for sg in ('k', 'w', 'gr'):
                    output[g][sg][i,j] = peaks[g][sg]
          
     
    return output


def make_gr(file):
    output = gen_gr()  
    with h5py.File(file, 'w') as f:
      
        for g in ('vblab', 'nb', 'vbprf', 'vcprf'):    
             f[g] = output[g]
        
        for g in ('F', 'A', 'C', 'E'):
            f.require_group(g)
            for sg in ('k', 'w', 'gr'):
                f[g][sg] = output[g][sg]





if __name__ == '__main__':
     
    """
    vb = .92
    k, wr, wi =  solveDispRel(nb=0.16, vb=vb, qb=5, mb=3, krange=(-10,10), guessFcn=bestGuess)
   
    fig, ax = plt.subplots( figsize = [4,4])
    
    vbline = k*vb
    
    wi = medfilt(wi, 25)
    
    #Filtering bc branch in solution makes an annoying bump
    wr = savgol_filter(wr, 501, 3)
    
    
    
    wifactor = 10
    ax.set_ylim((-1,10))
    ax.plot(k, wr + k*vb, '-', k, wi*wifactor, '--')
    ax.plot(k, vbline, '--')
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlim(-4,4)
    plt.show()
    
    peaks = classify_peaks(k, wr, wi)
    
    """
    file = os.path.join("C:", os.sep, "Users","Peter", "Desktop", "gr_save_C+4_He+1.hdf5")
    #file = '/Users/peter/Desktop/new_gr_save.hdf5'
    print(file)
    make_gr(file)
    
