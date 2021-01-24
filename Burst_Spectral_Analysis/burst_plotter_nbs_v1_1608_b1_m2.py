# this code will plot the results of time resolved spectral analysis of X-ray bursts :

import glob
from astropy.io import ascii
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Enter here the name of the source as in the burster_v3f.dat : 
source_name = '4U_1608-522'
# enter here the burstid of the burst you would like to fit :
bid = '1'
# Enter here the fit_method for which the results will be plotted :

fit_method = '3'

# maximum time range for the plot
max_time = 30.0
mine=input('Enter the minimum energy of the fits :')
mines = ":"+mine
maxe=input('Enter the maximum energy of the fits :')
maxes = ":"+maxe

burst_folder = '/home/hea/ownCloud/burst_characterization_v4/'
#burst_folder ='/home/hea/ownCloud/burst_characterization_v4/'
sfolder = '/home/hea/ownCloud/burst_characterization_v4/scripts/Burst_Spectral_Analysis/'
#sfolder = '/home/hea/ownCloud/burst_characterization_v3/scripts'

data = pd.read_csv(burst_folder+source_name+'/burst'+bid+'/'+source_name+'_sp_res_'+bid+fit_method+'_sbb.csv')


Source_Name_l=data['Source_Name_l']
BID_l=data['BID_l']
SID_l=data['SID_l']
OBS_ID_l=data['OBS_ID_l']
MJD_OBS_l=data['MJD_OBS_l']
exp_l=data['exp_l']
#b_time = MJD_OBS_l-MJD_OBS_l[0]*86400+exp_l/2.0
b_time=(MJD_OBS_l-MJD_OBS_l.min())*86400+exp_l/2.0-0.5
DATE_OBS_l=data['DATE_OBS_l']
exp_err = [exp_l/2.0,exp_l/2.0]
NH_l=data['NH_l']
BB_kT_l=data['BB_kT_l']
BBkT_err=[np.abs(data['min_BBkT_l']),data['max_BBkt_l']]
BB_Norm_l=data['BB_Norm_l']
min_BBNorm_l=['min_BBNorm_l']
max_BBNorm_l=['max_BBNorm_l']
min_BBkT_l=data['min_BBkT_l']
max_BBkt_l=data['max_BBkt_l']
BBNorm_err=[np.abs(data['min_BBNorm_l']),data['max_BBNorm_l']]
BB_flux_l=data['BB_flux_l']
BB_flux_tugba = BB_flux_l/1e-9
min_BBflux_l=data['min_BBflux_l']
max_BBflux_l=data['max_BBflux_l']
BB_flux_err = [data['min_BBflux_l']/1e-9,data['max_BBflux_l']/1e-9]
Bol_BBF_l=data['Bol_BBF_l']
Boll_tugba=Bol_BBF_l/1e-9
min_BolBBF_l=data['min_BolBBF_l']
max_BolBBF_l=data['max_BolBBF_l']

Bol_BBF_err = [data['min_BolBBF_l']/1e-9,data['max_BolBBF_l']/1e-9]
BB2_kT_l=data['BB2_kT_l']
BB2kT_err = [np.abs(data['min_2BBkT_l']),data['max_2BBkt_l']]
BB2_Norm_l=data['BB2_Norm_l']
BB2Norm_err=[np.abs(data['min_2BBNorm_l']),data['max_2BBNorm_l']]
BB2_flux_l=data['BB2_flux_l']
BB2_flux_err = [data['min_2BBflux_l']/1e-9,data['max_2BBflux_l']/1e-9]
Bol_2BBF_l=data['Bol_2BBF_l']
Bol_2BBF_err = [data['min_2BolBBF_l']/1e-9,data['max_2BolBBF_l']/1e-9]
fa_l=data['fa_l']
fa_err = [np.abs(data['min_fa_l']),data['max_fa_l']]
dbb_kT_l=data['dbb_kT_l']
dbbkT_err = [np.abs(data['min_dbbkT_l']),data['max_dbbkT_l']]
dbb_Norm_l=data['dbb_Norm_l']
dbb_Norm_err = [np.abs(data['min_dbbNorm_l']),data['max_dbbNorm_l']]
dbb_flux_l=data['dbb_flux_l']
dbbflux_err = [data['min_dbbflux_l']/1e-9,data['max_dbbflux_l']/1e-9]
sBB_kT_l=data['sBB_kT_l']
sBBkT_err = [np.abs(data['min_sBBkT_l']),data['max_sBBkt_l']]
sBB_Norm_l=data['sBB_Norm_l']
sBBNorm_err = [np.abs(data['min_sBBNorm_l']),data['max_sBBNorm_l']]
sBB_flux_l=data['sBB_flux_l']
sBBflux_err = [data['min_sBBflux_l']/1e-9,data['max_sBBflux_l']/1e-9]
RStat_l=data['RStat_l']
dof_l=data['dof_l']
rchi = RStat_l/dof_l


msize= 4.0
eline =1.5
#clf()
plt.clf()
#ioff()
plt.ioff()

if fit_method == '1':
    print('1=fixed background just BB free \n')
    print('Creating a four panel figure with Flux, kT, Norm, chi2')

    # Just to determine the plot ranges :
    sel_burst_time_range = np.where((b_time > 1.0) & (b_time < max_time) & (rchi > 0.0) & (BB_Norm_l > 0.0) & (BB_kT_l < 10.0) & (BB_kT_l > 0.05))

    print(sel_burst_time_range)
    kt_range_sel = BB_kT_l[sel_burst_time_range[0]]
    norm_range_sel = BB_Norm_l[sel_burst_time_range[0]]
    chi_range_sel = rchi[sel_burst_time_range[0]]
    min_range_kt = min(kt_range_sel)
    max_range_kt = max(kt_range_sel)
    min_range_norm = min(norm_range_sel)
    max_range_norm = max(norm_range_sel)
    max_range_chi = max(chi_range_sel)

    f, axarr = plt.subplots(4, sharex=True)
    plt.subplots_adjust(hspace=0.001)

    axarr[0].errorbar(b_time, Boll_tugba, xerr=exp_err, yerr=Bol_BBF_err, color='Black',linestyle='None', marker='.', label = 'BB Flux',ms=msize, elinewidth=eline)
    axarr[0].set_ylabel(r'Flux ($\times10^{-9}$)')
    #axarr[0].set_title(r' '+source_name+'  Burst:'+bid)
    axarr[0].set_ylim(-7.0,120.0)
    #axarr[0].plot([18.523,18.523],[0.01,160],linestyle='--',color='Red',lw=2) #--for red line
    axarr[0].grid(True)
    axarr[0].set_yticks([0.0,40.0,80.0,120.0])
    #axarr[0].set_title(r' '+'4U 1606-52  Burst:'+bid)


    axarr[1].errorbar(b_time, BB_kT_l, yerr=BBkT_err,xerr=exp_err,linestyle='None', marker='.', label = 'BB kT',ms=msize, elinewidth=eline,color='Black',)
    axarr[1].set_ylim(1.0,2.5)
    axarr[1].set_yticks([0.5,1.5,2.0])
    #axarr[1].set_ylim(min_range_kt-0.1*min_range_kt,max_range_kt+0.1*max_range_kt)
    #axarr[1].plot([18.523,18.523],[0.0,150],linestyle='--',color='Red',lw=2)
    axarr[1].grid(True)
    axarr[1].set_ylabel(r'kT (keV)')

    #axarr[2].set_ylim(min_range_norm-0.1*min_range_norm,max_range_norm+0.1*max_range_norm)
    axarr[2].set_ylim(-40.0,750.0)
    axarr[2].set_yticks([50.0,250.0,400.0])
    #axarr[2].set_ylim(min_range_norm-0.6*min_range_norm,max_range_norm+0.6*max_range_norm)
    #axarr[2].set_yscale('log')
    axarr[2].errorbar(b_time, BB_Norm_l, xerr=exp_err, yerr=BBNorm_err, linestyle='None',marker='.', label = 'BB Norm',ms=msize, elinewidth=eline,color='Black',)
    #axarr[2].plot([18.523,18.523],[0.01,1350],linestyle='--',color='Red',lw=2)
    axarr[2].grid(True)
    axarr[2].set_ylabel(r'Norm (R$^2$/D$_{10kpc})$' )
    
    axarr[3].errorbar(b_time, rchi, xerr=exp_err, linestyle='None', marker='.', label = 'BB Norm',ms=msize, elinewidth=eline,color='Black',)
    axarr[3].set_ylim(0.6,2.0)
    axarr[3].set_yticks([0.5,1.0,1.5])
    #axarr[3].set_ylim(0.5,max_range_chi+0.05)
    axarr[3].set_ylabel('$\chi^2/dof$')
    axarr[3].plot([0.0,1000.0],[1.0,1.0], linestyle='dashed', color='Green')
    #axarr[3].plot([18.523,18.523],[0.01,5],linestyle='--',color='Red',lw=2)
    axarr[3].grid(True)
    axarr[3].set_xlabel('Time (s)')    
    plt.xlim(1.0,max_time)
    
    
    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'.png',orientation='landscape')
    if not os.path.isdir(burst_folder+source_name+'/plots/'):
        print('No Plots folder found, creating a new one')
        os.makedirs(burst_folder+source_name+'/plots/')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'.png',orientation='landscape')

if fit_method == '3':
    print('3=fixed background one free BB and fa \n') 
    print('Creating a four panel figure with Flux, kT, Norm, fa, chi2')

    # Just to determine the plot ranges :
    sel_burst_time_range = np.where((np.abs(b_time) > 1.0) & (b_time < max_time) & (rchi > 0.0) & (BB_Norm_l > 0.0) & (BB_kT_l < 10.0) & (BB_kT_l > 0.05))

    kt_range_sel = BB_kT_l[sel_burst_time_range[0]]
    norm_range_sel = BB_Norm_l[sel_burst_time_range[0]]
    chi_range_sel = rchi[sel_burst_time_range[0]]
    min_range_kt = min(kt_range_sel)
    max_range_kt = max(kt_range_sel)
    min_range_norm = min(norm_range_sel)
    max_range_norm = max(norm_range_sel)
    max_range_chi = max(chi_range_sel)

    f, axarr = plt.subplots(5, sharex=True)  
    plt.subplots_adjust(hspace=0.001)
    
    axarr[0].errorbar(b_time, Bol_BBF_l/1e-9, xerr=exp_err, yerr=Bol_BBF_err,linestyle='None', marker='.', label = 'BB Flux',ms=msize, elinewidth=eline,color='Black',)
    axarr[0].set_ylabel(r'Flux ($\times10^{-9}$)')
    axarr[0].set_ylim(-10.0,240.0)
    #axarr[0].plot([18.523,18.523],[0.01,150],linestyle='--',color='Red',lw=2)
    axarr[0].set_title(r' '+'4U 1606-522 Burst:'+bid)

    axarr[1].errorbar(b_time, BB_kT_l, yerr=BBkT_err,xerr=exp_err,linestyle='None', marker='.', label = 'BB kT',ms=msize, elinewidth=eline,color='Black',)
    axarr[1].set_ylim(0.6,4.9)
    #axarr[1].plot([18.523,18.523],[0.0,150],linestyle='--',color='Red',lw=2)
    axarr[1].set_ylabel(r'kT (keV)')

    axarr[2].errorbar(b_time, BB_Norm_l, xerr=exp_err, yerr=BBNorm_err, linestyle='None',marker='.', label = 'BB Norm',ms=msize, elinewidth=eline,color='Black',)
    axarr[2].set_ylim(-100.0,900.0)
    #axarr[2].set_ylim(min_range_norm-0.1*min_range_norm,max_range_norm+0.1*max_range_norm)
    #axarr[2].set_yscale('log')
    #axarr[2].plot([18.523,18.523],[0.01,1150],linestyle='--',color='Red',lw=2)
    axarr[2].set_ylabel(r'Norm (R$^2$/D$_{10kpc})$' )

    axarr[3].errorbar(b_time, fa_l, xerr=exp_err, yerr=fa_err, linestyle='None',marker='.', label = 'Fa',ms=msize, elinewidth=eline,color='Black',)
    #axarr[3].set_yscale('log')
    axarr[3].set_ylim(0.5,2.4)
    axarr[3].plot([-6.0,max_time],[1.0,1.0],linestyle='--')
    #axarr[3].plot([18.523,18.523],[0.01,5],linestyle='--',color='Red',lw=2)
    axarr[3].set_ylabel(r'$f_{a}$' )

    axarr[4].errorbar(b_time, rchi, xerr=exp_err, linestyle='None', marker='.', label = 'BB Norm',ms=msize, elinewidth=eline,color='Black',)
    axarr[4].set_ylim(0.5, 1.7)
    axarr[4].set_ylabel('$\chi^2/dof$')
    axarr[4].plot([-0.1,1000.0],[1.0,1.0], linestyle='dashed', color='Green')
    #axarr[4].plot([18.523,18.523],[0.01,5],linestyle='--',color='Red',lw=2 )
    axarr[4].set_xlabel('Time (s)')    
    plt.xlim(1.8,max_time)
    
    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'_new_.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'_new_.png',orientation='landscape')
    if not os.path.isdir(burst_folder+source_name+'/plots/'):
        print('No Plots folder found, creating a new one')
        os.makedirs(burst_folder+source_name+'/plots/')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'_new_.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_'+mine+'_'+maxe+'_new_.png',orientation='landscape')



if (fit_method == '2') or (fit_method == '4'):
    if (fit_method == '2'):
        print('2=thawed background and one free BB \n')
    if (fit_method  == '4'):
        print('4=fixed background two BB \n') 

    print('Creating a four panel figure with Flux, kT, Norm, chi2 together with right axes showing the thawed parameter')

    # Just to determine the plot ranges :
    sel_burst_time_range = np.where((b_time > 1.0) & (b_time < max_time) & (rchi > 0.0) & (BB_Norm_l > 0.0) & (BB_kT_l < 10.0) & (BB_kT_l > 0.05))

    kt_range_sel = BB_kT_l[sel_burst_time_range[0]]
    norm_range_sel = BB_Norm_l[sel_burst_time_range[0]]
    chi_range_sel = rchi[sel_burst_time_range[0]]
    min_range_kt = min(kt_range_sel)
    max_range_kt = max(kt_range_sel)
    min_range_norm = min(norm_range_sel)
    max_range_norm = max(norm_range_sel)
    max_range_chi = max(chi_range_sel)

    # selection for the second BB

    bb2kt_range_sel=BB2_kT_l[sel_burst_time_range[0]]
    bb2norm_range_sel=BB2_Norm_l[sel_burst_time_range[0]]
    min_range_2kt = min(bb2kt_range_sel)
    max_range_2kt = max(bb2kt_range_sel)
    min_range_2norm = min(bb2norm_range_sel)
    max_range_2norm = max(bb2norm_range_sel)

    f, axarr = plt.subplots(4, sharex=True)
    plt.subplots_adjust(hspace=0.001)

    axarr[0].errorbar(b_time, Bol_BBF_l/1e-9, xerr=exp_err, yerr=Bol_BBF_err, linestyle='None', marker='.', label = 'BB Flux',ms=msize, elinewidth=eline, color='Black')
    axarr[0].set_ylabel(r'Flux ($\times10^{-9}$)')
    axarr[0].set_title(r' '+source_name+'  Burst #'+bid)
    axarr[0].set_ylim(-1.0,90.0)
    ax0 = axarr[0].twinx()
    if fit_method == '2':
        ax0.errorbar(b_time, Bol_2BBF_l/1e-9, yerr=Bol_2BBF_err,xerr=exp_err,linestyle='None', marker='.', label = 'BB 2 Flux',color='Red',ms=msize, elinewidth=eline)
        ax0.set_ylabel(r'Flux$_{\rm{Cool}}$', color='Red')
        ax0.tick_params(axis='y', colors='red')
        ax0.spines['right'].set_color('red')
    if fit_method == '4':
        ax0.errorbar(b_time, dbb_flux_l/1e-9, yerr=dbbflux_err,xerr=exp_err,linestyle='None', marker='.', label = 'DBB Flux',color='green',ms=msize, elinewidth=eline)
        ax0.errorbar(b_time, sBB_flux_l/1e-9, yerr=sBBflux_err,xerr=exp_err,linestyle='None', marker='.', label = 'sBB Flux',color='green',ms=msize, elinewidth=eline)
        ax0.set_ylabel(r'Flux', color='green')
        ax0.tick_params(axis='y', colors='green') 
        ax0.spines['right'].set_color('green')
    
    axarr[1].errorbar(b_time, BB_kT_l, yerr=BBkT_err,xerr=exp_err,linestyle='None', marker='.', label = 'BB kT',ms=msize, elinewidth=eline)
    axarr[1].set_ylim(-1.0,3.0)
    axarr[1].set_ylabel(r'kT (keV)')
    ax1 = axarr[1].twinx()
    if fit_method == '2':
        ax1.errorbar(b_time, BB2_kT_l, yerr=BB2kT_err,xerr=exp_err,linestyle='None', marker='.', label = 'BB 2 kT',color='Red',ms=msize, elinewidth=eline)
        ax1.set_ylabel(r'kT$_{\rm{Cool}}$', color='Red')
        ax1.tick_params(axis='y', colors='red')
        ax1.spines['right'].set_color('red')
    if fit_method == '4':
        ax1.errorbar(b_time, dbb_kT_l, yerr=dbb_kT_l,xerr=exp_err,linestyle='None', marker='.', label = 'DBB Flux',color='green',ms=msize, elinewidth=eline)
        ax1.errorbar(b_time, sBB_kT_l, yerr=sBB_kT_l,xerr=exp_err,linestyle='None', marker='.', label = 'sBB Flux',color='green',ms=msize, elinewidth=eline)
        ax1.set_ylabel(r'kT$_{\rm{Cool}}$', color='green')
        ax1.tick_params(axis='y', colors='green')
        ax1.spines['right'].set_color('green')

    axarr[2].set_ylim(min_range_norm-0.1*min_range_norm,max_range_norm+0.1*max_range_norm)
    axarr[2].set_ylim(-100.0,1000.0)
    #axarr[2].set_yscale('log')
    axarr[2].errorbar(b_time, BB_Norm_l, xerr=exp_err, yerr=BBNorm_err, linestyle='None',marker='.', label = 'BB Norm',ms=msize, elinewidth=eline)
    axarr[2].set_ylabel(r'Norm (R$^2$/D$_{10kpc})$' )
    ax2 = axarr[2].twinx()
    if fit_method == '2':
        ax2.errorbar(b_time, BB2_Norm_l, yerr=BB2_Norm_l,xerr=exp_err,linestyle='None', marker='.', label = 'BB 2 Norm',color='Red',ms=msize, elinewidth=eline)
        ax2.set_ylabel(r'Norm$_{\rm{Cool}}$', color='Red')
        ax2.tick_params(axis='y', colors='red')
        ax2.spines['right'].set_color('red')
    if fit_method == '4':
        ax2.errorbar(b_time, dbb_Norm_l, yerr=dbb_Norm_err,xerr=exp_err,linestyle='None', marker='.', label = 'DBB Flux',color='darkgrey',ms=msize, elinewidth=eline)
        ax2.errorbar(b_time, sBB_Norm_l, yerr=sBBNorm_err,xerr=exp_err,linestyle='None', marker='.', label = 'sBB Flux',color='green',ms=msize, elinewidth=eline)
        ax2.set_ylabel(r'Norm$_{\rm{Cool}}$', color='green')
        ax2.tick_params(axis='y', colors='green')
        ax2.spines['right'].set_color('green')
    
    axarr[3].errorbar(b_time, rchi, xerr=exp_err, linestyle='None', marker='.', label = 'BB Norm',ms=msize, elinewidth=eline)
    axarr[3].set_ylim(0.5,max_range_chi+max_range_chi*0.1)
    axarr[3].set_ylabel('$\chi^2/dof$')
    axarr[3].plot([0.0,1000.0],[1.0,1.0], linestyle='dashed', color='Green')
    axarr[3].set_xlabel('Time (s)')    
    plt.xlim(0.5,max_time)

    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'_new_.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/burst'+bid+'/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'_new_.png',orientation='landscape')
    if not os.path.isdir(burst_folder+source_name+'/plots/'):
        print('No Plots folder found, creating a new one')
        os.makedirs(burst_folder+source_name+'/plots/')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'_new_.pdf',orientation='landscape')
    plt.savefig(burst_folder+source_name+'/plots/figure'+source_name+'_b'+bid+'_f'+fit_method+'_new_'+mine+'_'+maxe+'_new_.png',orientation='landscape')











