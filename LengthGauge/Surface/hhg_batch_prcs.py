# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 12:48:37 2019

@author: BaykushevaDenitsaRan
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 20:11:27 2019

@author: denitsab
"""

import numpy as np 
from numpy.fft import fft , fftfreq, rfft, rfftfreq, hfft, fftshift
import matplotlib.pyplot as plt
from os import listdir
from os.path import isdir
def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n
import re

def retrieve_numvals(linein):
    return np.array(re.findall(r"[-+]?\d*\.*\d+", linein)).astype(float)
freq_0 = 1.5498/27.211*800/7500.0

masterdir_ = './calcscan2/'
dirs_ = [d for d in listdir(masterdir_) if ( ("I_" in d ) and isdir(masterdir_+d) )]
calcID=''
readinfo = True
for d in dirs_:
    plt.close("all")

    dir_ = masterdir_+d+'/'
    
    fs = [f for f in listdir(dir_) if "interband" in f]
    if( len(fs) >0 ):
        if( not (len(calcID) == 0) ):
            fname_inter = dir_+'interband_dipole_'+calcID+'.dat'
            fname_intra = dir_+'intraband_dipole_'+calcID+'.dat'
            fname_ncc   = dir_+'cpopulation_'+calcID+'.dat'
            fname_pvc   = dir_+'cvcoherence_'+calcID+'.dat'
            fname_info  = dir_+'info_'+calcID+'.dat'
        else:
            fname_inter = dir_ + [f for f in listdir(dir_) if "interband"  in f][0]
            fname_intra = dir_ + [f for f in listdir(dir_) if "intraband"  in f][0]
            fname_ncc   = dir_ + [f for f in listdir(dir_) if "population" in f][0]
            fname_pvc   = dir_ + [f for f in listdir(dir_) if "coherence"  in f][0]
            fname_info  = dir_ + [f for f in listdir(dir_) if (("info"       in f) and (".dat" in f))][0]
        if(readinfo):
            fid = open(fname_info, "U")
            fid.seek(0,0)
            bzcut = 'XX'
            for iline, line in enumerate(fid):
                if ("Time spent in integration routine" in line):
                    tcpu = retrieve_numvals(line)[0]
                    line1=line
                    print(retrieve_numvals(line))
                    print(tcpu)
                elif( "cycles" in line):
                    ncyc = retrieve_numvals(line)[0]
                elif("dephasing time" in line):
                    T2f = retrieve_numvals(line)[0]
                elif("kx-dimension" in line):
                    nkx= retrieve_numvals(line)[0]
                elif("ky-dimension" in line):
                    nky =  retrieve_numvals(line)[0]
                elif("t-domain" in line):
                    ntpts = retrieve_numvals(line)[0]
                elif("center wavelength" in line):
                    lambda_cen= retrieve_numvals(line)[0]
                elif("W/cm2" in line):
                    ITWcm2 = retrieve_numvals(line)[0]
                elif("f0" in line):
                    f0dk = retrieve_numvals(line)[1]
                elif("sigma" in line):
                    sigp = retrieve_numvals(line)[0]
                elif("hmin" in line):
                    hmin = retrieve_numvals(line)[0]
                elif("hmax" in line):
                    hmax= retrieve_numvals(line)[0]
                elif("h0" in line):
                    h0 = retrieve_numvals(line)[1]
                elif("Brillouin" in line):
                    bzcut = retrieve_numvals(line)[0]
            fid.close()
            calcinfo = "I_"+str(ITWcm2)+"TWcm2_ncyc_"+str(ncyc)+"_wav_"+str(lambda_cen)+"nm_T2_"+\
            str(T2f)+"per_nkx_"+str(nkx)+"_nky_"+str(nky)+"_f0_"+str(f0dk)+"_sig_"+\
            str(sigp)+"_zvode_"+str(hmin)+"_"+str(hmax)+"_"+str(h0)+'_BZ_'+str(bzcut)
        else:
            calcinfo=""
        
        
        ncc = np.loadtxt(fname_ncc, skiprows=0)
        ncc = ncc[:,1]
        pvc = np.loadtxt(fname_pvc, skiprows=0)
        pvcR = pvc[:,1]
        pvcI = pvc[:,2]
        dip_array = np.loadtxt(fname_inter,skiprows=0)
        tvec = dip_array[:,0]
        nt=len(tvec)
        dt = tvec[1]-tvec[0]
        dx_inter = dip_array[:,1]#/dt/dt
        dy_inter = dip_array[:,3]#/dt/dt
        
        dip_array = np.loadtxt(fname_intra,skiprows=0)
        dx_intra = dip_array[:,1]
        dy_intra = dip_array[:,3]
        
        dx_intra0 = dx_intra
        dx_inter0 = dx_inter
        dy_inter0 = dy_inter
        dy_intra0 = dy_intra
        
        # In[]
        
        filterfun = np.hanning(nt)
        #filterfun = 1.0
        subdc = 1.0
        dx_intra = (dx_intra - subdc*np.mean(dx_intra))*filterfun
        dy_intra = (dy_intra - subdc*np.mean(dy_intra))*filterfun
        dx_inter = (dx_inter - subdc*np.mean(dx_inter))*filterfun
        dy_inter = (dy_inter - subdc*np.mean(dy_inter))*filterfun
        
        
        
        NFFT = 16*nextpow2(len(tvec))
        NFFT = len(tvec)
        
        
        Yx_intra = fft(dx_intra, n = NFFT)[0:NFFT/2]/nt
        Yy_intra = fft(dy_intra, n = NFFT)[0:NFFT/2]/nt
        Yx_inter = fft(dx_inter, n = NFFT)[0:NFFT/2]/nt
        Yy_inter = fft(dy_inter, n = NFFT)[0:NFFT/2]/nt
        
        f_max = 0.5*1.0/( tvec[1] - tvec[0] )
        
        freq = (fftfreq(d=tvec[1]-tvec[0], n = NFFT)*2.*np.pi)
        #freq  = np.linspace(0.0, f_max, NFFT/2)*2.0*np.pi
        n_freq =freq/freq_0
        
        n_freq = n_freq[0:NFFT/2]
        
        # In[]
        xmin = 0
        xmax = 40
        plt.figure(facecolor='white')
        plt.semilogy(n_freq,  (np.abs(Yx_intra)**2))
        plt.semilogy(n_freq,  (np.abs(Yx_inter)**2))
        plt.semilogy(n_freq,  (np.abs(Yx_intra+Yx_inter)**2))
        plt.legend(['intraband','interband','coherent sum'])
        plt.xlabel(r'$\mathrm{frequency} / \omega_0$' )
        plt.ylabel(r'$\mathrm{intensity} / \mathrm{a.u.}$' )
        plt.title('parallel emission')
        plt.xlim([xmin,xmax])
        plt.tight_layout()
        plt.savefig(masterdir_+calcinfo+"_parallel.png")
        
        plt.figure(facecolor='white')
        plt.semilogy(n_freq,  (np.abs(Yy_intra)**2))
        plt.semilogy(n_freq,  (np.abs(Yy_inter)**2))
        plt.semilogy(n_freq,  (np.abs(Yy_intra+Yy_inter)**2))
        plt.legend(['intraband','interband','coherent sum'])
        plt.xlabel(r'$\mathrm{frequency} / \omega_0$' )
        plt.ylabel(r'$\mathrm{intensity} / \mathrm{a.u.}$' )
        plt.title('perpendicular emission')
        plt.xlim([xmin,xmax])
        plt.tight_layout()
        plt.savefig(masterdir_+calcinfo+"_orthogonal.png")
        
        plt.figure(facecolor='white')
        plt.semilogy(n_freq,  (np.abs(Yx_intra+Yx_inter)**2), 'b')
        plt.semilogy(n_freq,  (np.abs(Yy_intra+Yy_inter)**2), 'r')
        plt.legend(['parallel', 'perpendicular'])
        plt.xlabel(r'$\mathrm{frequency} / \omega_0$' )
        plt.ylabel(r'$\mathrm{intensity} / \mathrm{a.u.}$' )
        plt.title('total emission')
        plt.xlim([xmin,xmax])
        plt.tight_layout()
        plt.savefig(masterdir_+calcinfo+"_total.png")
        
        #raise_window()
        
#        # In[]
#        fig, ax = plt.subplots(facecolor='white')
#        ax2 = ax.twinx()
#        
#        ax2.plot(tvec, ncc)
#        ax2.legend("population")
#        ax2.set_ylabel(r"$\pi_{vc} (t)$ / arb.u.")
#        
#        ax.plot(tvec, pvcR,"r")
#        ax.plot( tvec,pvcI,"g")
#        
#        ax.set_xlabel(r"$t$/au")
#        ax.set_ylabel(r"$n_{cc}(t)$/ arb. u.")
#        #ax.set_ylim([pvcR.min(),pvcR.max()])
#        ax.legend(["real part coherence","imag. part coherence"])
#        fig.tight_layout()
#        
#        # In[]
#        fig, ax2 = plt.subplots(facecolor='white')
#        
#        ax2.plot(tvec, dx_inter0,"r")
#        ax2.plot( tvec,dy_inter0,"c")
#        
#        ax2.set_xlabel(r"$t$/au")
#        ax2.set_ylabel(r"$J_{er}(t)$/ arb. u.")
#        #ax.set_ylim([pvcR.min(),pvcR.max()])
#        ax2.legend(["interband X","interband Y"])
#        
#        fig, ax = plt.subplots(facecolor='white')
#        ax.plot(tvec, dx_intra0,"r")
#        ax.plot( tvec,dy_intra0,"g")
#        
#        ax.set_xlabel(r"$t$/au")
#        ax.set_ylabel(r"$J_{ra}(t)$/ arb. u.")
#        #ax.set_ylim([pvcR.min(),pvcR.max()])
#        ax.legend(["intraband X","intraband Y"])
#        fig.tight_layout()