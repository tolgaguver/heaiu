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

# Create the dataframe to be filled later :

Source_Name_l=[]
BID_l=[]
SID_l=[]
OBS_ID_l=[]
MJD_OBS_l=[]
DATE_OBS_l=[]
exp_l=[]
NH_l=[]
BB_kT_l=[]
min_BBkT_l=[]
max_BBkt_l=[]
BB_Norm_l=[]
min_BBNorm_l=[]
max_BBNorm_l=[]
BB_flux_l=[]
min_BBflux_l=[]
max_BBflux_l=[]
Bol_BBF_l=[]
min_BolBBF_l=[]
max_BolBBF_l=[]
BB2_kT_l=[]
min_2BBkT_l=[]
max_2BBkt_l=[]
BB2_Norm_l=[]
min_2BBNorm_l=[]
max_2BBNorm_l=[]
BB2_flux_l=[]
min_2BBflux_l=[]
max_2BBflux_l=[]
Bol_2BBF_l=[]
min_2BolBBF_l=[]
max_2BolBBF_l=[]
fa_l=[]
min_fa_l=[]
max_fa_l=[]
dbb_kT_l=[]
min_dbbkT_l=[]
max_dbbkT_l=[]
dbb_Norm_l=[]
min_dbbNorm_l=[]
max_dbbNorm_l=[]
dbb_flux_l=[]
min_dbbflux_l=[]
max_dbbflux_l=[]
sBB_kT_l=[]
min_sBBkT_l=[]
max_sBBkt_l=[]
sBB_Norm_l=[]
min_sBBNorm_l=[]
max_sBBNorm_l=[]
sBB_flux_l=[]
min_sBBflux_l=[]
max_sBBflux_l=[]
RStat_l=[]
dof_l=[]

# Enter here the name of the source as in the burster_v3f.dat : 
source_name = '4U_1608-522'
# enter here the burstid of the burst you would like to fit :
bid = input('enter here the burstid of the burst you would like to fit : ')
#bid = '4'
mine=input('Enter the minimum energy of the fits :')
mines = ":"+mine
maxe=input('Enter the maximum energy of the fits :')
maxes = maxe+":"

folder = '/home/hea/ownCloud/burst_characterization_v3/'
sfolder = '/home/hea/ownCloud/burst_characterization_v3/scripts/'

#folder = '/Users/tolga/ownCloud/burst_characterization_v3/'
#sfolder = '/Users/tolga/ownCloud/burst_characterization_v3/scripts/'


# Read the persistent state analysis (pre burst tbabs*DISKBB+BBODYRAD fit)
pers_file = folder+source_name+'/pers_results.dat'
pers_data  = pd.read_csv(pers_file)
snh = pers_data['NH']
diskbb_temp = pers_data['disk_kT']
diskbb_norm = pers_data['disk_norm']
sbb_kt=pers_data['bb_kT']
sbb_norm=pers_data['bb_norm']
pers_rstat = pers_data['chi']
pers_dof = pers_data['dof']

sel_pre_burst = np.where(pers_data['BID']==int(bid))
ssnh = snh[sel_pre_burst[0]]
sdiskbb_temp = diskbb_temp[sel_pre_burst[0]]
sdiskbb_norm = diskbb_norm[sel_pre_burst[0]]
ssbb_kt = sbb_kt[sel_pre_burst[0]]
ssbb_norm = sbb_norm[sel_pre_burst[0]]
spers_rstat= pers_rstat[sel_pre_burst[0]]
spers_dof = pers_dof[sel_pre_burst[0]]
print('These are the values I get from persistent state analysis :')
print('NH ='+str(ssnh.values[0]))
print('Disk BB kT = '+str(sdiskbb_temp.values[0]))
print('Disk BB Norm = '+str(sdiskbb_norm.values[0]))
print('Surface BB kT = '+str(ssbb_kt.values[0]))
print('Surface BB Norm = '+str(ssbb_norm.values[0]))
print('Reduced Chi2 of = '+str(spers_rstat.values[0]/spers_dof.values[0]))

print('available fit_methods \n')
print('4=fixed background two BB \n') 

#fit_method = input('Please enter your fit preference : ')
fit_method = '4'
burst_folder=folder+source_name+'/burst'+bid+'/'
bkg_folder = folder+source_name+'/burst'+bid+'/pers_analysis/'
bkgfile = glob.glob(bkg_folder+'*3c50*.pha')
sp_list = np.array(glob.glob(burst_folder+'c*.pha'))
pha_count = len(sp_list)

set_stat("chi2xspecvar")
set_covar_opt("sigma",1.0)
set_conf_opt('numcores', 10)
set_conf_opt("max_rstat",250.0)
set_covar_opt('sigma',1.0)


for i in range(len(sp_list)-1):
    sp_final = sp_list[i]
    print(sp_final)
    sp_hdu = fits.open(str(sp_final))
    print('read src spec: '+str(sp_final))
    mjdobs = sp_hdu[1].header['MJD-OBS']
    date_obsi = sp_hdu[1].header['DATE-OBS']
    exposure = sp_hdu[1].header['EXPOSURE']
    obsid = sp_hdu[1].header['OBS_ID']
    sid = sp_final.split('/')[7].split('_')[0][1:]
    date_obs = str(Time(date_obsi,format='isot', scale='utc'))
    print(date_obs)
    print(obsid)
    object = sp_hdu[1].header['OBJECT']
    Source_Name_l.append(object)
    BID_l.append(bid)
    SID_l.append(sid)
    OBS_ID_l.append(obsid)
    MJD_OBS_l.append(mjdobs)
    DATE_OBS_l.append(date_obs)
    exp_l.append(exposure)
    NH_l.append(ssnh.values[0])

    if exposure == 0.0 :
        print('this spectrum does not have any exposure')
        BB_kT_l.append(0)
        min_BBkT_l.append(0)
        max_BBkt_l.append(0)
        BB_Norm_l.append(0)
        min_BBNorm_l.append(0)
        max_BBNorm_l.append(0)
        BB_flux_l.append(0)
        min_BBflux_l.append(0)
        max_BBflux_l.append(0)
        Bol_BBF_l.append(0)
        min_BolBBF_l.append(0)
        max_BolBBF_l.append(0)
        BB2_kT_l.append(0)
        min_2BBkT_l.append(0)
        max_2BBkt_l.append(0)
        BB2_Norm_l.append(0)
        min_2BBNorm_l.append(0)
        max_2BBNorm_l.append(0)
        BB2_flux_l.append(0)
        min_2BBflux_l.append(0)
        max_2BBflux_l.append(0)
        Bol_2BBF_l.append(0)
        min_2BolBBF_l.append(0)
        max_2BolBBF_l.append(0)
        fa_l.append(0)
        min_fa_l.append(0)
        max_fa_l.append(0)
        dbb_kT_l.append(0)
        min_dbbkT_l.append(0)
        max_dbbkT_l.append(0)
        dbb_Norm_l.append(0)
        min_dbbNorm_l.append(0)
        max_dbbNorm_l.append(0)
        dbb_flux_l.append(0)
        min_dbbflux_l.append(0)
        max_dbbflux_l.append(0)
        sBB_kT_l.append(0)
        min_sBBkT_l.append(0)
        max_sBBkt_l.append(0)
        sBB_Norm_l.append(0)
        min_sBBNorm_l.append(0)
        max_sBBNorm_l.append(0)
        sBB_flux_l.append(0)
        min_sBBflux_l.append(0)
        max_sBBflux_l.append(0)
        RStat_l.append(0)
        dof_l.append(0)
        print('Finished:') 
        print(sp_final)
        continue
    
    load_pha(1, str(sp_final),use_errors=True)
    load_arf(1, sfolder+'ni_xrcall_onaxis_v1.02.arf')
    load_rmf(1, sfolder+'nicer_v1.02.rmf')
    print('Ignoring : '+mines+' '+maxes)
    ignore(mines+','+maxes)
    print('This script only subtracts ni3c50 background')
    load_bkg(1, bkgfile[0])
    subtract()
    print('Grouping the data to have at least 50 counts per channel')
    group_counts(1, 50)

    # first let's do a global source definition :
    set_source(xstbabs.tb*(xsdiskbb.dbb+xsbbodyrad.sbb))
    tb.nH=ssnh.values[0]
    dbb.Tin = sdiskbb_temp.values[0]
    dbb.norm = sdiskbb_norm.values[0]
    sbb.kT = ssbb_kt.values[0]
    sbb.norm = ssbb_norm.values[0]
    freeze(tb.nH)
    freeze(dbb.Tin)
    freeze(dbb.norm)
    freeze(sbb.kT)
    freeze(sbb.norm)
    #fit()
    initial_rstat = sum(calc_chisqr())/len(calc_chisqr()) 
    if (initial_rstat < (spers_rstat.values[0])/spers_dof.values[0]):
        print('Current chi2 : '+str(initial_rstat))
        print('Persistent chi2 : '+str(spers_rstat.values[0]/spers_dof.values[0]))
        print('Deviation from the persistent emission is small I will thaw the parameters and refit the data to save the best fit values:')
        thaw(dbb.Tin)
        thaw(dbb.norm)
        thaw(sbb.kT)
        thaw(sbb.norm)
        fit()
        if get_fit_results().dof <= 0.0:
            print('The degree of freedom is very small we need to skip this spectrum:')
            dbb_kT_l.append(0)
            min_dbbkT_l.append(0)
            max_dbbkT_l.append(0)
            dbb_Norm_l.append(0)
            min_dbbNorm_l.append(0)
            max_dbbNorm_l.append(0)
            dbb_flux_l.append(0)
            min_dbbflux_l.append(0)
            max_dbbflux_l.append(0)
            sBB_kT_l.append(0)
            min_sBBkT_l.append(0)
            max_sBBkt_l.append(0)
            sBB_Norm_l.append(0)
            min_sBBNorm_l.append(0)
            max_sBBNorm_l.append(0)
            sBB_flux_l.append(0)
            min_sBBflux_l.append(0)
            max_sBBflux_l.append(0)
            BB_kT_l.append(0)
            min_BBkT_l.append(0)
            max_BBkt_l.append(0)
            BB_Norm_l.append(0)
            min_BBNorm_l.append(0)
            max_BBNorm_l.append(0)
            BB_flux_l.append(0)
            min_BBflux_l.append(0)
            max_BBflux_l.append(0)
            Bol_BBF_l.append(0)
            min_BolBBF_l.append(0)
            max_BolBBF_l.append(0)
            BB2_kT_l.append(0)
            min_2BBkT_l.append(0)
            max_2BBkt_l.append(0)
            BB2_Norm_l.append(0)
            min_2BBNorm_l.append(0)
            max_2BBNorm_l.append(0)
            BB2_flux_l.append(0)
            min_2BBflux_l.append(0)
            max_2BBflux_l.append(0)
            Bol_2BBF_l.append(0)
            min_2BolBBF_l.append(0)
            max_2BolBBF_l.append(0)
            fa_l.append(0)
            min_fa_l.append(0)
            max_fa_l.append(0)
            RStat_l.append(0)
            dof_l.append(get_fit_results().dof)
            print('Finished:') 
            print(sp_final)
            continue
        covar()
        chi = get_fit_results().statval
        dof = get_fit_results().dof
        parvals = np.array(get_covar_results().parvals)
        parnames = np.array(get_covar_results().parnames)
        parmins = np.array(get_covar_results().parmins)
        parmaxes = np.array(get_covar_results().parmaxes)
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)                
                
        if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:

            print('It seems like you have unconstrained parameters in the continuum model therefore we cant calculate the errors')
            print('No flux will be reported')                        
            dbb_flux_l.append(0)
            max_dbbflux_l.append(0)
            min_dbbflux_l.append(0)
            sBB_flux_l.append(0)
            max_sBBflux_l.append(0)
            min_sBBflux_l.append(0)
                          
            print('Parameter Errors will be written as 0')
            dbb_kT_l.append(get_fit_results().parvals[0])
            min_dbbkT_l.append(0)
            max_dbbkT_l.append(0)
            dbb_Norm_l.append(get_fit_results().parvals[1])
            min_dbbNorm_l.append(0)
            max_dbbNorm_l.append(0)
            sBB_kT_l.append(get_fit_results().parvals[2])
            min_sBBkT_l.append(0)
            max_sBBkt_l.append(0)
            sBB_Norm_l.append(get_fit_results().parvals[3])
            min_sBBNorm_l.append(0)
            max_sBBNorm_l.append(0)

            #elif ((None in cparmins) == False and (None in cparmaxes) == False) or ((0 in cparmaxes) == False and (0 in cparmins) == False):
        else:
            print('The parameters are constrained well calculating errors')
            #matrix = get_covar_results().extra_output
            #is_all_zero = np. all((matrix > 0))
            #if is_all_zero: 
            sample2=sample_flux(dbb,float(mine),float(maxe), num=100, correlated=False,confidence=68)
            dbb_flux_l.append(sample2[1][0]) 
            max_dbbflux_l.append(sample2[1][1]-sample2[1][0])
            min_dbbflux_l.append(sample2[1][0]-sample2[1][2])
            sample3=sample_flux(sbb,float(mine),float(maxe), num=100, correlated=False,confidence=68)
            sBB_flux_l.append(sample3[1][0]) 
            max_sBBflux_l.append(sample3[1][1]-sample3[1][0])
            min_sBBflux_l.append(sample3[1][0]-sample3[1][2]) 
            # Parameter errors will be written as they are :
            dbb_kT_l.append(get_covar_results().parvals[0])
            min_dbbkT_l.append(get_covar_results().parmins[0])
            max_dbbkT_l.append(get_covar_results().parmaxes[0])
            dbb_Norm_l.append(get_covar_results().parvals[1])
            min_dbbNorm_l.append(get_covar_results().parmins[1])
            max_dbbNorm_l.append(get_covar_results().parmaxes[1])
            sBB_kT_l.append(get_covar_results().parvals[2])
            min_sBBkT_l.append(get_covar_results().parmins[2])
            max_sBBkt_l.append(get_covar_results().parmaxes[2])
            sBB_Norm_l.append(get_covar_results().parvals[3])
            min_sBBNorm_l.append(get_covar_results().parmins[3])
            max_sBBNorm_l.append(get_covar_results().parmaxes[3])
        BB_kT_l.append(0)
        min_BBkT_l.append(0)
        max_BBkt_l.append(0)
        BB_Norm_l.append(0)
        min_BBNorm_l.append(0)
        max_BBNorm_l.append(0)
        BB_flux_l.append(0)
        min_BBflux_l.append(0)
        max_BBflux_l.append(0)
        Bol_BBF_l.append(0)
        min_BolBBF_l.append(0)
        max_BolBBF_l.append(0)
        BB2_kT_l.append(0)
        min_2BBkT_l.append(0)
        max_2BBkt_l.append(0)
        BB2_Norm_l.append(0)
        min_2BBNorm_l.append(0)
        max_2BBNorm_l.append(0)
        BB2_flux_l.append(0) 
        max_2BBflux_l.append(0)
        min_2BBflux_l.append(0)
        Bol_2BBF_l.append(0)
        min_2BolBBF_l.append(0)
        max_2BolBBF_l.append(0)
        fa_l.append(0)
        min_fa_l.append(0)
        max_fa_l.append(0)
	
        RStat_l.append(chi)
        dof_l.append(dof)

        plot_data()
        x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.05
        x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.05)
        ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
        ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)

        plot_fit_delchi(1,clearwindow=True, color='Black')
        fig=plt.gcf()
        ax1,ax2=fig.axes
        ax1.set_title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax2.set_xscale('log')
        ax2.set_xlabel('Energy [keV]', fontsize=14)
        ax1.set_ylabel('Counts/sec/keV', fontsize=14)
        ax2.set_ylabel('Sigma', fontsize=14)
        ax1.set_xlim(x_min,x_max)
        ax2.set_xlim(x_min,x_max)
        plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
        plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')
        plot_fit(1,clearwindow=True,xlog=True,ylog=True, color='Black')
        plot_model_component("tb*dbb", replot=False, overplot=True, color='Green')
        plot_model_component("tb*sbb", replot=False, overplot=True, color='Red')
        plt.title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)
        plt.xlabel('Energy [keV]', fontsize=14)
        plt.ylabel('Counts/sec/keV', fontsize=14)
        plt.xlim(x_min,x_max)
        plt.ylim(ymin,ymax)
        plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
        plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')
        print('We are done with just fitting the continuum, Skipping to the next spectrum')
        continue
    else:
        print('Current chi2 : '+str(initial_rstat))
        print('Persistent chi2 : '+str(spers_rstat.values[0]/spers_dof.values[0]))
        print('Simple Persistent Emission Does Not Fit the data we need to add components')
	
        print('This is 4=fixed background two BB') 

        set_source(xstbabs.tb*(xsbbodyrad.bb+xsdiskbb.dbb+xsbbodyrad.sbb))
        tb.nH=ssnh.values[0]
        dbb.Tin = sdiskbb_temp.values[0]
        dbb.norm = sdiskbb_norm.values[0]
        sbb.kT = ssbb_kt.values[0]
        sbb.norm = ssbb_norm.values[0]
        freeze(tb.nH)
        freeze(dbb.Tin)
        freeze(dbb.norm)
        freeze(sbb.kT)
        freeze(sbb.norm)
        dbb_kT_l.append(diskbb_temp.values[0])
        min_dbbkT_l.append(0)
        max_dbbkT_l.append(0)
        dbb_Norm_l.append(diskbb_norm.values[0])
        min_dbbNorm_l.append(0)
        max_dbbNorm_l.append(0)
        dbb_flux_l.append(0)
        min_dbbflux_l.append(0)
        max_dbbflux_l.append(0)
        sBB_kT_l.append(ssbb_kt.values[0])
        min_sBBkT_l.append(0)
        max_sBBkt_l.append(0)
        sBB_Norm_l.append(ssbb_norm.values[0])
        min_sBBNorm_l.append(0)
        max_sBBNorm_l.append(0)
        sBB_flux_l.append(0)
        min_sBBflux_l.append(0)
        max_sBBflux_l.append(0)
        fa_l.append(0)
        min_fa_l.append(0)
        max_fa_l.append(0)
        bb.kt=0.5
        set_xsabund('wilm')
        bb.norm = 180.3
        set_method("moncar")
        fit()
        set_method("levmar")
        fit()
        chi = get_fit_results().statval
        dof = get_fit_results().dof
        if get_fit_results().dof <= 0.0:
            print('The degree of freedom is very small we need to skip this spectrum:')
            BB_kT_l.append(0)
            min_BBkT_l.append(0)
            max_BBkt_l.append(0)
            BB_Norm_l.append(0)
            min_BBNorm_l.append(0)
            max_BBNorm_l.append(0)
            BB_flux_l.append(0)
            min_BBflux_l.append(0)
            max_BBflux_l.append(0)
            Bol_BBF_l.append(0)
            min_BolBBF_l.append(0)
            max_BolBBF_l.append(0)
            BB2_kT_l.append(0)
            min_2BBkT_l.append(0)
            max_2BBkt_l.append(0)
            BB2_Norm_l.append(0)
            min_2BBNorm_l.append(0)
            max_2BBNorm_l.append(0)
            Bol_2BBF_l.append(0)
            max_2BolBBF_l.append(0)
            min_2BolBBF_l.append(0)
            BB2_flux_l.append(0) 
            max_2BBflux_l.append(0)
            min_2BBflux_l.append(0)
            RStat_l.append(0)
            dof_l.append(get_fit_results().dof)
            print('Finished:') 
            print(sp_final)
            continue
        if (get_fit_results().rstat>=1.3):
            print('Fit is not acceptable trying to add a second blackbody')
            set_source(tb*(bb+xsbbodyrad.bb2+dbb+sbb))
            tb.nH=ssnh.values[0]
            dbb.Tin = sdiskbb_temp.values[0]
            dbb.norm = sdiskbb_norm.values[0]
            sbb.kT = ssbb_kt.values[0]
            sbb.norm = ssbb_norm.values[0]
            freeze(tb.nH)
            freeze(dbb.Tin)
            freeze(dbb.norm)
            freeze(sbb.kT)
            freeze(sbb.norm)
            bb2.kT=0.2
            bb2.norm=1000
            set_method("moncar")
            fit()
            set_method("levmar")
            covar()
            chi = get_fit_results().statval
            dof = get_fit_results().dof
            RStat_l.append(chi)
            dof_l.append(dof)
            if get_fit_results().dof <= 0.0:
                print('The degree of freedom is very small we need to skip this spectrum:')
                BB_kT_l.append(0)
                min_BBkT_l.append(0)
                max_BBkt_l.append(0)
                BB_Norm_l.append(0)
                min_BBNorm_l.append(0)
                max_BBNorm_l.append(0)
                BB_flux_l.append(0)
                min_BBflux_l.append(0)
                max_BBflux_l.append(0)
                Bol_BBF_l.append(0)
                min_BolBBF_l.append(0)
                max_BolBBF_l.append(0)
                BB2_kT_l.append(0)
                min_2BBkT_l.append(0)
                max_2BBkt_l.append(0)
                BB2_Norm_l.append(0)
                min_2BBNorm_l.append(0)
                max_2BBNorm_l.append(0)
                Bol_2BBF_l.append(0)
                max_2BolBBF_l.append(0)
                min_2BolBBF_l.append(0)
                BB2_flux_l.append(0) 
                max_2BBflux_l.append(0)
                min_2BBflux_l.append(0)
                print('Finished:') 
                print(sp_final)
                continue
            parvals = np.array(get_covar_results().parvals)
            parnames = np.array(get_covar_results().parnames)
            parmins = np.array(get_covar_results().parmins)
            parmaxes = np.array(get_covar_results().parmaxes)
                     
            # now the model unabsorbed fluxes :
            covar()
            cparmins = np.array(get_covar_results().parmins)
            cparmaxes = np.array(get_covar_results().parmaxes)
            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
                print('It seems like you have unconstrained parameters we cant calculate errors')
                print('No flux will be reported')
                BB_flux_l.append(0)
                max_BBflux_l.append(0)
                min_BBflux_l.append(0)
                Bol_BBF_l.append(0)
                min_BolBBF_l.append(0)
                max_BolBBF_l.append(0)
                # Parameter Errors will be written as 0:
                BB_kT_l.append(get_covar_results().parvals[0])
                min_BBkT_l.append(0)
                max_BBkt_l.append(0)
                BB_Norm_l.append(get_covar_results().parvals[1])
                min_BBNorm_l.append(0)
                max_BBNorm_l.append(0)
                BB2_kT_l.append(get_covar_results().parvals[2])
                min_2BBkT_l.append(0)
                max_2BBkt_l.append(0)
                BB2_Norm_l.append(get_covar_results().parvals[3])
                min_2BBNorm_l.append(0)
                max_2BBNorm_l.append(0)
                Bol_2BBF_l.append(0)
                max_2BolBBF_l.append(0)
                min_2BolBBF_l.append(0)
                BB2_flux_l.append(0) 
                max_2BBflux_l.append(0)
                min_2BBflux_l.append(0)
                    
            else:
                #elif ((None in cparmins) == False and (None in cparmaxes) == False) or ((0 in cparmaxes) == False and (0 in cparmins) == False):
                print('The parameters are constrained well calculating errors')
                # matrix = get_covar_results().extra_output
                # is_all_zero = np. all((matrix > 0))
                # if is_all_zero: 
                sample1=sample_flux(bb,float(mine),float(maxe), num=100, correlated=False,confidence=68)
                BB_flux_l.append(sample1[1][0]) 
                max_BBflux_l.append(sample1[1][1]-sample1[1][0])
                min_BBflux_l.append(sample1[1][0]-sample1[1][2])
                sample2=sample_flux(bb2,float(mine),float(maxe), num=100, correlated=False,confidence=68)
                BB2_flux_l.append(sample2[1][0]) 
                max_2BBflux_l.append(sample2[1][1]-sample2[1][0])
                min_2BBflux_l.append(sample2[1][0]-sample2[1][2])
                # Parameter errors will be written as they are 
                BB_kT_l.append(get_covar_results().parvals[0])
                min_BBkT_l.append(get_covar_results().parmins[0])
                max_BBkt_l.append(get_covar_results().parmaxes[0])
                BB_Norm_l.append(get_covar_results().parvals[1])
                min_BBNorm_l.append(get_covar_results().parmins[1])
                max_BBNorm_l.append(get_covar_results().parmaxes[1])
                BB2_kT_l.append(get_covar_results().parvals[2])
                min_2BBkT_l.append(get_covar_results().parmins[2])
                max_2BBkt_l.append(get_covar_results().parmaxes[2])
                BB2_Norm_l.append(get_covar_results().parvals[3])
                min_2BBNorm_l.append(get_covar_results().parmins[3])
                max_2BBNorm_l.append(get_covar_results().parmaxes[3])
                           
                # Now the Bolometric Fluxes :
                sample_bol=sample_flux(bb,0.01,200.0, num=100, correlated=False,confidence=68)
                Bol_BBF_l.append(sample_bol[1][0]) 
                max_BolBBF_l.append(sample_bol[1][1]-sample_bol[1][0])
                min_BolBBF_l.append(sample_bol[1][0]-sample_bol[1][2])
                # Now the Bolometric Fluxes :
                sample_bol2=sample_flux(bb2,0.01,200.0, num=100, correlated=False,confidence=68)
                Bol_2BBF_l.append(sample_bol2[1][0]) 
                max_2BolBBF_l.append(sample_bol2[1][1]-sample_bol2[1][0])
                min_2BolBBF_l.append(sample_bol2[1][0]-sample_bol2[1][2])
       
            plot_data()
            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.1
            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.1)
            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)

            plot_fit_delchi(1,clearwindow=True, color='Black')
            fig=plt.gcf()
            ax1,ax2=fig.axes
            ax1.set_title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            ax2.set_xscale('log')
            ax2.set_xlabel('Energy [keV]', fontsize=14)
            ax1.set_ylabel('Counts/sec/keV', fontsize=14)
            ax2.set_ylabel('Sigma', fontsize=14)
            ax1.set_xlim(x_min,x_max)
            ax2.set_xlim(x_min,x_max)
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
                
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')
            
            plot_fit(1,clearwindow=True,xlog=True,ylog=True, color='Black')
            
            plot_model_component("tb*dbb", replot=False, overplot=True, color='Green')                        
            plot_model_component("tb*sbb", replot=False, overplot=True, color='Red')                        
            plot_model_component("tb*bb", replot=False, overplot=True, color='Black')
            plot_model_component("tb*bb2", replot=False, overplot=True, color='Blue')       
            plt.title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)                        
            plt.xlabel('Energy [keV]', fontsize=14)                        
            plt.ylabel('Counts/sec/keV', fontsize=14)                        
            plt.xlim(x_min,x_max)                        
            plt.ylim(ymin,ymax)
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')
            continue
        if (get_fit_results().rstat<1.5):
            print('Current chi2 : '+str(get_fit_results().rstat))
            print('It looks like just one blackbody is enough so I will just write that result')
            BB2_kT_l.append(0)
            min_2BBkT_l.append(0)
            max_2BBkt_l.append(0)
            BB2_Norm_l.append(0)
            min_2BBNorm_l.append(0)
            max_2BBNorm_l.append(0)
            Bol_2BBF_l.append(0)
            max_2BolBBF_l.append(0)
            min_2BolBBF_l.append(0)
            BB2_flux_l.append(0) 
            max_2BBflux_l.append(0)
            min_2BBflux_l.append(0)
            chi = get_fit_results().statval
            dof = get_fit_results().dof
            RStat_l.append(chi)
            dof_l.append(get_fit_results().dof)
            # now the model unabsorbed fluxes :
            covar()
            cparmins = np.array(get_covar_results().parmins)
            cparmaxes = np.array(get_covar_results().parmaxes)
            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:

                print('It seems like you have unconstrained parameters we cant calculate errors')
                print('No flux will be reported')
                BB_flux_l.append(0)
                max_BBflux_l.append(0)
                min_BBflux_l.append(0)
                Bol_BBF_l.append(0)
                min_BolBBF_l.append(0)
                max_BolBBF_l.append(0)
                # Parameter Errors will be written as 0:
                BB_kT_l.append(get_covar_results().parvals[0])
                min_BBkT_l.append(0)
                max_BBkt_l.append(0)
                BB_Norm_l.append(get_covar_results().parvals[1])
                min_BBNorm_l.append(0)
                max_BBNorm_l.append(0)
                            
            else:
                #elif ((None in cparmins) == False and (None in cparmaxes) == False) or ((0 in cparmaxes) == False and (0 in cparmins) == False):
                print('The parameters are constrained well calculating errors')
                
                sample1=sample_flux(bb,float(mine),float(maxe), num=100, correlated=False,confidence=68)
                BB_flux_l.append(sample1[1][0]) 
                max_BBflux_l.append(sample1[1][1]-sample1[1][0])
                min_BBflux_l.append(sample1[1][0]-sample1[1][2])
                print('Parameter errors will be written as they are')
                BB_kT_l.append(get_covar_results().parvals[0])
                min_BBkT_l.append(get_covar_results().parmins[0])
                max_BBkt_l.append(get_covar_results().parmaxes[0])
                BB_Norm_l.append(get_covar_results().parvals[1])
                min_BBNorm_l.append(get_covar_results().parmins[1])
                max_BBNorm_l.append(get_covar_results().parmaxes[1])
                           
                # Now the Bolometric Fluxes :
                sample_bol=sample_flux(bb,0.01,200.0, num=100, correlated=False,confidence=68)
                Bol_BBF_l.append(sample_bol[1][0]) 
                max_BolBBF_l.append(sample_bol[1][1]-sample_bol[1][0])
                min_BolBBF_l.append(sample_bol[1][0]-sample_bol[1][2])
                                  
                        
            plot_data()
            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.1
            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.1)
            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)

            plot_fit_delchi(1,clearwindow=True, color='Black')
            fig=plt.gcf()
            ax1,ax2=fig.axes
            ax1.set_title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            ax2.set_xscale('log')
            ax2.set_xlabel('Energy [keV]', fontsize=14)
            ax1.set_ylabel('Counts/sec/keV', fontsize=14)
            ax2.set_ylabel('Sigma', fontsize=14)
            ax1.set_xlim(x_min,x_max)
            ax2.set_xlim(x_min,x_max)
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
                
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')
            
            plot_fit(1,clearwindow=True,xlog=True,ylog=True, color='Black')
            
            plot_model_component("tb*dbb", replot=False, overplot=True, color='Green')                        
            plot_model_component("tb*sbb", replot=False, overplot=True, color='Red')                        
            plot_model_component("tb*bb", replot=False, overplot=True, color='Black')                       
            plt.title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)                        
            plt.xlabel('Energy [keV]', fontsize=14)                        
            plt.ylabel('Counts/sec/keV', fontsize=14)                        
            plt.xlim(x_min,x_max)                        
            plt.ylim(ymin,ymax)
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.pdf',orientation='landscape', papertype='a4')
            plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_'+mine+'_'+maxe+'.png',orientation='landscape', papertype='a4')

items = {'Source_Name' : Source_Name_l, 'BID' : BID_l, 'SID': SID_l,'OBS_ID' : OBS_ID_l,'MJD-OBS' : MJD_OBS_l,'DATE-OBS' : OBS_ID_l,'exp' :exp_l,'NH' : NH_l,
'BB_kT' : BB_kT_l,'min_BBkT' : min_BBkT_l,'max_BBkt' : max_BBkt_l,'BB_Norm' : BB_Norm_l,'min_BBNorm' : min_BBNorm_l,'max_BBNorm' : max_BBNorm_l,
'BB_flux' : BB_flux_l,'min_BBflux' : min_BBflux_l,
'max_BBflux' : max_BBflux_l,'Bol_BBF' : Bol_BBF_l,'min_BolBBF' : min_BolBBF_l,'max_BolBBF' : max_BolBBF_l,'2BB_kT' : BB2_kT_l,
'min_2BBkT' : min_2BBkT_l,'max_2BBkt' : max_2BBkt_l,'2BB_Norm' : BB2_Norm_l,
'min_2BBNorm' : min_2BBNorm_l,'max_2BBNorm' : max_2BBNorm_l,'2BB_flux' : BB2_flux_l,'min_2BBflux' : min_2BBflux_l,'max_2BBflux' : max_2BBflux_l,'Bol_2BBF' : Bol_2BBF_l,'min_2BolBBF' : min_2BolBBF_l,
'max_2BolBBF' : max_2BolBBF_l,'fa' : fa_l,'min_fa' : min_fa_l,'max_fa' : max_fa_l,'dbb_kT' : dbb_kT_l,'min_dbbkT' : min_dbbkT_l,'max_dbbkT' : max_dbbkT_l,'dbb_Norm' : dbb_Norm_l,'min_dbbNorm' : min_dbbNorm_l,
'max_dbbNorm' : max_dbbNorm_l,'dbb_flux' : dbb_flux_l,'min_dbbflux' : min_dbbflux_l,'max_dbbflux' : max_dbbflux_l,
'sBB_kT' : sBB_kT_l,'min_sBBkT' : min_sBBkT_l,'max_sBBkt' : max_sBBkt_l,'sBB_Norm' : sBB_Norm_l,'min_sBBNorm' : min_sBBNorm_l,
'max_sBBNorm' : max_sBBNorm_l,'sBB_flux' : sBB_flux_l,
'min_sBBflux' : min_sBBflux_l,'max_sBBflux' : max_sBBflux_l,
'RStat' : RStat_l,'dof' : dof_l}

#print(len(Source_Name_l))
print('Sourcename:',len(Source_Name_l))
print('BID:',len(BID_l))
print('SID:',len(SID_l))
print('OBSID:',len(OBS_ID_l))
print('MJD:',len(MJD_OBS_l))
print('exp:',len(exp_l))
print('NH:',len(NH_l))
print('BB_kT:',len(BB_kT_l))
print('min_BBkT:',len(min_BBkT_l))
print('max_BBkt:',len(max_BBkt_l))
print('BB_Norm:',len(BB_Norm_l))
print('min_BBNorm:',len(min_BBNorm_l))
print('max_BBNorm:',len(max_BBNorm_l))
print('BB_flux:',len(BB_flux_l))
print('min_BBflux:',len(min_BBflux_l))
print('max_BBflux:',len(max_BBflux_l))
print('Bol_BBF:',len(Bol_BBF_l))
print('min_BolBBF:',len(min_BolBBF_l))
print('max_BolBBF:',len(max_BolBBF_l))
print('BB2_kT:',len(BB2_kT_l))
print('min_2BBkT:',len(min_2BBkT_l))
print('max_2BBkt:',len(max_2BBkt_l))
print('BB2_Norm:',len(BB2_Norm_l))
print('min_2BBNorm:',len(min_2BBNorm_l))
print('max_2BBNorm:',len(max_2BBNorm_l))
print('BB2_flux:',len(BB2_flux_l))
print('min_2BBflux:',len(min_2BBflux_l))
print('max_2BBflux:',len(max_2BBflux_l))
print('Bol_2BBF:',len(Bol_2BBF_l))
print('min_2BolBBF:',len(min_2BolBBF_l))
print('max_2BolBBF:',len(max_2BolBBF_l))
print('fa:',len(fa_l))
print('min_fa:',len(min_fa_l))
print('max_fa:',len(max_fa_l))
print('dbb_kT:',len(dbb_kT_l))
print('min_dbbkT:',len(min_dbbkT_l))
print('max_dbbkT:',len(max_dbbkT_l))
print('dbb_Norm:',len(dbb_Norm_l))
print('min_dbbNorm:',len(min_dbbNorm_l))
print('max_dbbNorm:',len(max_dbbNorm_l))
print('dbb_flux:',len(dbb_flux_l))
print('min_dbbflux:',len(min_dbbflux_l))
print('max_dbbflux:',len(max_dbbflux_l))
print('sBB_kT:',len(sBB_kT_l))
print('min_sBBkT:',len(min_sBBkT_l))                             
print('max_sBBkt:',len(max_sBBkt_l))
print('sBB_Norm:',len(sBB_Norm_l))
print('min_sBBNorm:',len(min_sBBNorm_l))
print('max_sBBNorm:',len(max_sBBNorm_l))
print('sBB_flux:',len(sBB_flux_l))
print('min_sBBflux:',len(min_sBBflux_l))
print('max_sBBflux:',len(max_sBBflux_l))
print('Rstat:',len(RStat_l))
print('Dof:',len(dof_l))
#for i in check_list:
#print(i,str(len(i)))

df=pd.DataFrame.from_dict(items)
df.transpose()
#df.to_excel(burst_folder+source_name+'_sp_res_'+bid+'_'+fit_method+'.xlsx') 
#df.to_csv(burst_folder+source_name+'_sp_res_'+bid+fit_method+'.csv') 

sorted_df = df.sort_values(by='MJD-OBS')
sorted_df.to_excel(burst_folder+source_name+'_sp_res_'+bid+'_'+fit_method+'_'+mine+'_'+maxe+'.xlsx') 
sorted_df.to_csv(burst_folder+source_name+'_sp_res_'+bid+fit_method+'_'+mine+'_'+maxe+'.csv') 
                

		

