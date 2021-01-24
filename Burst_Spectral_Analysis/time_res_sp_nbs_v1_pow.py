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
col_names = ["Source_Name_l", "BID_l", "SID_l", "OBS_ID_l", "MJD_OBS_l",\
             "DATE_OBS_l", "exp_l", "NH_l", "BB_kT_l", "min_BBkT_l","max_BBkt_l",\
             "BB_Norm_l", "min_BBNorm_l", "max_BBNorm_l", "BB_flux_l", \
             "min_BBflux_l", "max_BBflux_l", "Bol_BBF_l", "min_BolBBF_l", \
             "max_BolBBF_l", "BB2_kT_l", "min_2BBkT_l", "max_2BBkt_l", "BB2_Norm_l",\
             "min_2BBNorm_l", "max_2BBNorm_l", "BB2_flux_l", "min_2BBflux_l",\
             "max_2BBflux_l", "Bol_2BBF_l", "min_2BolBBF_l", "max_2BolBBF_l",\
             "fa_l", "min_fa_l", "max_fa_l", "dbb_kT_l", "min_dbbkT_l", "max_dbbkT_l",\
             "dbb_Norm_l", "min_dbbNorm_l", "max_dbbNorm_l", "dbb_flux_l", "min_dbbflux_l",\
             "max_dbbflux_l", "sBB_kT_l", "min_sBBkT_l", "max_sBBkt_l", "sBB_Norm_l",\
             "min_sBBNorm_l", "max_sBBNorm_l", "sBB_flux_l", "min_sBBflux_l", "max_sBBflux_l",\
             "RStat_l", "dof_l", "mine", "maxe", "r_mine", "r_maxe"]

def print_all_lengths(data_frame):
    for key in data_frame:
        print("%s: %d" %(key, len(data_frame[key])))
    return

def create_plots(components, x_max, x_min, ymax, ymin):
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
    plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_pow.pdf',orientation='landscape', papertype='a4')
    plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_pow.png',orientation='landscape', papertype='a4')
    plot_fit(1,clearwindow=True,xlog=True,ylog=True, color='Black')
    if "dbb" in components:
        plot_model_component("tb*dbb", replot=False, overplot=True, color='Green')
    if "sbb" in components:
        plot_model_component("tb*sbb", replot=False, overplot=True, color='Red')
    if "bb" in components:
        plot_model_component("tb*bb", replot=False, overplot=True, color='Black')
    if "bb2" in components:
        plot_model_component("tb*bb2", replot=False, overplot=True, color='Blue')
    if "fa_dbb" in components:        
        plot_model_component("tb*fa*dbb", replot=False, overplot=True, color='Green')
    if "fa_sbb" in components:
        plot_model_component("tb*fa*sbb", replot=False, overplot=True, color='Red')                 
    plt.title(source_name+' BID:'+bid+' MJD:'+str(mjdobs)+' SID:'+sid)
    plt.xlabel('Energy [keV]', fontsize=14)
    plt.ylabel('Counts/sec/keV', fontsize=14)
    plt.xlim(x_min,x_max)
    plt.ylim(ymin,ymax)
    plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_pow.pdf',orientation='landscape', papertype='a4')
    plt.savefig(burst_folder+sid+'_b'+bid+'_'+fit_method+'_full_comp_pow.png',orientation='landscape', papertype='a4')
    print('Skipping to the next spectrum')
    

# Enter here the name of the source as in the burster_v3f.dat : 
source_name = '4U_1608-522'
# enter here the burstid of the burst you would like to fit :
bid = input('Enter the BID of the burst you would like to analyze :' ) #'1'
print('This is the burst that you are fitting :'+str(bid))
mine=input('Enter the minimum energy of the fits :')
mines = ":"+mine
maxe=input('Enter the maximum energy of the fits :')
maxes = maxe+":"


folder = '/home/hea/ownCloud/burst_characterization_v4/'
sfolder = '/home/hea/ownCloud/burst_characterization_v4/scripts/'
#folder = '/Users/tolga/ownCloud/burst_characterization_v3/'
#sfolder = '/Users/tolga/ownCloud/burst_characterization_v3/scripts/'

# Read the persistent state analysis (pre burst tbabs*DISKBB+POW fit)
pers_file = folder+source_name+'/pers_results.dat'
pers_data  = pd.read_csv(pers_file)
#X=np.asarray(pers_data).astype(np.float64)
#snh = pers_data['NH']
#diskbb_temp = pers_data['disk_kT']
#diskbb_norm = pers_data['disk_norm']
#sbb_kt=pers_data['bb_kT']
#sbb_norm=pers_data['bb_norm']
#pers_rstat = pers_data['chi']
#pers_dof = pers_data['dof']

sel_pre_burst = np.where(pers_data['BID']==int(bid))
ssnh = pers_data['NH'][sel_pre_burst[0]]
sdiskbb_temp = pers_data['disk_kT'][sel_pre_burst[0]]
sdiskbb_norm = pers_data['disk_norm'][sel_pre_burst[0]]
ssbb_kt = pers_data['bb_kT'][sel_pre_burst[0]]
ssbb_norm = pers_data['bb_norm'][sel_pre_burst[0]]
spers_rstat=  pers_data['chi'][sel_pre_burst[0]]
spers_dof = pers_data['dof'][sel_pre_burst[0]]

print('These are the values I get from persistent state analysis :')
print('NH ='+str(ssnh.values[0]))
print('Disk BB kT = '+str(sdiskbb_temp.values[0]))
print('Disk BB Norm = '+str(sdiskbb_norm.values[0]))
print('Surface BB kT = '+str(ssbb_kt.values[0]))
print('Surface BB Norm = '+str(ssbb_norm.values[0]))
print('Reduced Chi2 of = '+str(spers_rstat.values[0]/spers_dof.values[0]))

print('available fit_methods \n')
print('1=fixed background just BB free \n')
print('2=thawed background and one free BB \n')
print('3=fixed background one free BB and fa \n') 
print('4=fixed background two BB \n') 

fit_method = input('Please enter your fit preference : ')
# fit_method = '1'



burst_folder=folder+source_name+'/burst'+bid+'/'
bkg_folder = folder+source_name+'/burst'+bid+'/'
bkgfile = glob.glob(bkg_folder+'*3c50*.pha.pi')
sp_list = np.array(glob.glob(burst_folder+'c*.pha'))
pha_count = len(sp_list)

data_frame = {col: [0]*len(sp_list) for col in col_names}

set_stat("chi2xspecvar")
set_covar_opt("sigma",1.0)
set_conf_opt('numcores', 10)
set_conf_opt("max_rstat",250.0)
set_covar_opt('sigma',1.0)

# for i in range(len(sp_list)-1):
for i in range(len(sp_list)):
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
        data_frame["Source_Name_l"][i] = (object)
        data_frame["BID_l"][i] = (bid)
        data_frame["SID_l"][i] = (sid)
        data_frame["OBS_ID_l"][i] = (obsid)
        data_frame["MJD_OBS_l"][i] = (mjdobs)
        data_frame["DATE_OBS_l"][i] = (date_obs)
        data_frame["exp_l"][i] = (exposure)
        data_frame["NH_l"][i] = (ssnh.values[0])
        data_frame["mine"][i] = mine
        data_frame["maxe"][i] = maxe

        
        if exposure == 0.0 :
            print('this spectrum does not have any exposure')
            print('Finished:') 
            print(sp_final)
            continue
        # print(date_obsi)
        else:
            load_pha(1, str(sp_final),use_errors=True)
            load_arf(1, sfolder+'nicer-arf-consim135o-teamonly-array52.arf')
            load_rmf(1, sfolder+'nicer-rmf6s-teamonly-array52.rmf')
            load_bkg(1, bkgfile[0],use_errors=True)
            subtract()
            print('This script only subtracts ni3c50 background')
            print('Grouping the data to have at least 50 counts per channel')
            group_counts(1, 50)
            print('Ignoring : '+mines+' '+maxes)
            ignore(mines+','+maxes)
            if sum(get_counts(1,filter=True)) <= 100:
                print('This spectrum does not have enough counts to fit after background subtraction. Skipping')
                continue
            real_range = get_filter().split(":")
            print(real_range[0], real_range[1])
            data_frame["r_mine"][i] = real_range[0]
            data_frame["r_maxe"][i] = real_range[1]
   
    
            # first let's do a global source definition :
            set_source(xstbabs.tb*(xsdiskbb.dbb+xspowerlaw.sbb))
            tb.nH=ssnh.values[0]
            dbb.Tin = sdiskbb_temp.values[0]
            dbb.norm = sdiskbb_norm.values[0]
            sbb.PhoIndex = ssbb_kt.values[0]
            sbb.norm = ssbb_norm.values[0]
            freeze(tb.nH)
            freeze(dbb.Tin)
            freeze(dbb.norm)
            freeze(sbb.PhoIndex)
            freeze(sbb.norm)
            #fit()
            initial_rstat = sum(calc_chisqr())/len(calc_chisqr()) 
            if (initial_rstat < (spers_rstat.values[0]+3.50)/spers_dof.values[0]):
                    print('Current chi2 : '+str(initial_rstat))
                    print('Persistent chi2 : '+str(spers_rstat.values[0]/spers_dof.values[0]))
                    print('Deviation from the persistent emission is small I will thaw the parameters and refit the data to save the best fit values:')
                    thaw(dbb.Tin)
                    thaw(dbb.norm)
                    thaw(sbb.PhoIndex)
                    thaw(sbb.norm)
                    fit()
                    if get_fit_results().dof <= 0.0:
                        print('The degree of freedom is very small we need to skip this spectrum:')    
                        data_frame["dof_l"][i] = (get_fit_results().dof)
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
                    #if get_covar_results().parmins[1] == None: 
                    #        min_dbbNorm_l[i] = (0)
                    #if get_covar_results().parmins[1] != None: 
                    #        min_dbbNorm_l[i] = (get_covar_results().parmins[1])
                    
                    if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
    
                            print('It seems like you have unconstrained parameters we cant calculate errors')
                            print('No flux will be reported')                              
                            data_frame["dbb_kT_l"][i] = (get_covar_results().parvals[0])
                            data_frame["dbb_Norm_l"][i] = (get_covar_results().parvals[1])
                            data_frame["sBB_kT_l"][i] = (get_covar_results().parvals[2])
                            data_frame["sBB_Norm_l"][i] = (get_covar_results().parvals[3])
    
                    if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                            print('The parameters are constrained well calculating errors')
                            #matrix = get_covar_results().extra_output
                            #is_all_zero = np. all((matrix >= 0))
                            #if is_all_zero: 
                            sample2=sample_flux(dbb,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                            data_frame["dbb_flux_l"][i] = (sample2[1][0]) 
                            data_frame["max_dbbflux_l"][i] = (sample2[1][1]-sample2[1][0])
                            data_frame["min_dbbflux_l"][i] = (sample2[1][0]-sample2[1][2])
                            sample3=sample_flux(sbb,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                            data_frame["sBB_flux_l"][i] = (sample3[1][0]) 
                            data_frame["max_sBBflux_l"][i] = (sample3[1][1]-sample3[1][0])
                            data_frame["min_sBBflux_l"][i] = (sample3[1][0]-sample3[1][2]) 
                            # Parameter errors will be written as they are :
                            data_frame["dbb_kT_l"][i] = (get_covar_results().parvals[0])
                            data_frame["min_dbbkT_l"][i] = (get_covar_results().parmins[0])
                            data_frame["max_dbbkT_l"][i] = (get_covar_results().parmaxes[0])
                            data_frame["dbb_Norm_l"][i] = (get_covar_results().parvals[1])
                            data_frame["min_dbbNorm_l"][i] = (get_covar_results().parmins[1])
                            data_frame["max_dbbNorm_l"][i] = (get_covar_results().parmaxes[1])
                            data_frame["sBB_kT_l"][i] = (get_covar_results().parvals[2])
                            data_frame["min_sBBkT_l"][i] = (get_covar_results().parmins[2])
                            data_frame["max_sBBkt_l"][i] = (get_covar_results().parmaxes[2])
                            data_frame["sBB_Norm_l"][i] = (get_covar_results().parvals[3])
                            data_frame["min_sBBNorm_l"][i] = (get_covar_results().parmins[3])
                            data_frame["max_sBBNorm_l"][i] = (get_covar_results().parmaxes[3])
                    
                    
                    data_frame["RStat_l"][i] = (chi)
                    data_frame["dof_l"][i] = (dof)
    
                    plot_data()
                    x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.05
                    x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.05)
                    ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                    ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                    
                    create_plots(["dbb", "sbb"], x_max, x_min, ymax, ymin)
                    print('Skipping to the next spectrum')
                    continue
            else:
                    print('Current chi2 : '+str(initial_rstat))
                    print('Persistent chi2 : '+str(spers_rstat.values[0]/spers_dof.values[0]))
                    print('Simple Persistent Emission Does Not Fit the data we need to add components')
    
                    if (fit_method == '1') or (fit_method == '4'):
                            if (fit_method == '1'):
                                    print('This is 1=fixed background just BB free')
                            else:
                                    print('This is 4=fixed background two BB') 
                            set_source(xstbabs.tb*(xsbbodyrad.bb+xsdiskbb.dbb+xspowerlaw.sbb))
                            tb.nH=ssnh.values[0]
                            dbb.Tin = sdiskbb_temp.values[0]
                            dbb.norm = sdiskbb_norm.values[0]
                            sbb.PhoIndex = ssbb_kt.values[0]
                            sbb.norm = ssbb_norm.values[0]
                            freeze(tb.nH)
                            freeze(dbb.Tin)
                            freeze(dbb.norm)
                            freeze(sbb.PhoIndex)
                            freeze(sbb.norm)
                            
                            data_frame["dbb_kT_l"][i] = (sdiskbb_temp.values[0])
                            data_frame["dbb_Norm_l"][i] = (sdiskbb_norm.values[0])
                            data_frame["sBB_Norm_l"][i] = (ssbb_norm.values[0])
    
                                                    
                    elif (fit_method == '2'):
                            print('2=thawed background and one free BB')
                            set_source(xstbabs.tb*(xsbbodyrad.bb+xsdiskbb.dbb+xspowerlaw.sbb))
                            tb.nh=ssnh.values[0]
                            dbb.Tin = sdiskbb_temp.values[0]
                            dbb.norm = sdiskbb_norm.values[0]
                            sbb.PhoIndex = ssbb_kt.values[0]
                            sbb.norm = ssbb_norm.values[0]
                            thaw(dbb.Tin)
                            thaw(dbb.norm)
                            thaw(sbb.PhoIndex)
                            thaw(sbb.norm)
                            freeze(tb.nH)
                    elif fit_method == '3':
                            print('3=fixed background one free BB and fa') 
                            set_source(xstbabs.tb*(xsbbodyrad.bb+scale1d.fa*(xsdiskbb.dbb+xspowerlaw.sbb)))
                            tb.nH=ssnh.values[0]
                            dbb.Tin = sdiskbb_temp.values[0]
                            dbb.norm = sdiskbb_norm.values[0]
                            sbb.PhoIndex = ssbb_kt.values[0]
                            sbb.norm = ssbb_norm.values[0]
                            freeze(tb.nH)
                            freeze(dbb.Tin)
                            freeze(dbb.norm)
                            freeze(sbb.PhoIndex)
                            freeze(sbb.norm)
                            set_par (fa.c0, val=1.0, min=0.9)
                            
                            data_frame["dbb_kT_l"][i] = (sdiskbb_temp.values[0])
                            data_frame["dbb_Norm_l"][i] = (sdiskbb_norm.values[0])
                            data_frame["sBB_kT_l"][i] = (ssbb_kt.values[0])
                            data_frame["sBB_Norm_l"][i] = (ssbb_norm.values[0])
                                             
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
                            
                            data_frame["dof_l"][i] = (get_fit_results().dof)
                            print('Finished:') 
                            print(sp_final)
                            continue
                        
                    if (fit_method == '4') and (get_fit_results().rstat>=1.5):
                            print('Fit is not acceptable trying to add a second blackbody')
                            set_source(tb*(bb+xsbbodyrad.bb2+dbb+sbb))
                            tb.nH=ssnh.values[0]
                            dbb.Tin = sdiskbb_temp.values[0]
                            dbb.norm = sdiskbb_norm.values[0]
                            sbb.PhoIndex = ssbb_kt.values[0]
                            sbb.norm = ssbb_norm.values[0]
                            freeze(tb.nH)
                            freeze(dbb.Tin)
                            freeze(dbb.norm)
                            freeze(sbb.PhoIndex)
                            freeze(sbb.norm)
                            bb2.kT=0.2
                            fit()
                            covar()
                            chi = get_fit_results().statval
                            dof = get_fit_results().dof
                            parvals = np.array(get_covar_results().parvals)
                            parnames = np.array(get_covar_results().parnames)
                            parmins = np.array(get_covar_results().parmins)
                            parmaxes = np.array(get_covar_results().parmaxes)
                            #BB_kT_l[i] = (get_covar_results().parvals[0])
                            #min_BBkT_l[i] = (get_covar_results().parmins[0])
                            #max_BBkt_l[i] = (get_covar_results().parmaxes[0])
                            #BB_Norm_l[i] = (get_covar_results().parvals[1])
                            #min_BBNorm_l[i] = (get_covar_results().parmins[1])
                            #max_BBNorm_l[i] = (get_covar_results().parmaxes[1])
                            #BB2_kT_l[i] = (get_covar_results().parvals[2])
                            #min_2BBkT_l[i] = (get_covar_results().parmins[2])
                            #max_2BBkt_l[i] = (get_covar_results().parmaxes[2])
                            #BB2_Norm_l[i] = (get_covar_results().parvals[3])
                            #min_2BBNorm_l[i] = (get_covar_results().parmins[3])
                            #max_2BBNorm_l[i] = (get_covar_results().parmaxes[3])
                    
                            # now the model unabsorbed fluxes :
                            covar()
                            cparmins = np.array(get_covar_results().parmins)
                            cparmaxes = np.array(get_covar_results().parmaxes)
                            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
                                    print('It seems like you have unconstrained parameters we cant calculate errors')
                                    print('No flux will be reported')
                                    
                                    # Parameter Errors will be written as 0:
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[1])
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[2])
                                    data_frame["BB2_kT_l"][i] = (get_covar_results().parvals[3])
                                    data_frame["BB2_Norm_l"][i] = (get_covar_results().parvals[3])
                                    
                            
                            if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                                    print('The parameters are constrained well calculating errors')
                                    #matrix = get_covar_results().extra_output
                                    #is_all_zero = np. all((matrix >= 0))
                                    #if is_all_zero: 
                                    sample1=sample_flux(bb,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                                    data_frame["BB_flux_l"][i] = (sample1[1][0]) 
                                    data_frame["max_BBflux_l"][i] = (sample1[1][1]-sample1[1][0])
                                    data_frame["min_BBflux_l"][i] = (sample1[1][0]-sample1[1][2])
                                    sample2=sample_flux(bb2,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                                    data_frame["BB2_flux_l"][i] = (sample2[1][0]) 
                                    data_frame["max_2BBflux_l"][i] = (sample2[1][1]-sample2[1][0])
                                    data_frame["min_2BBflux_l"][i] = (sample2[1][0]-sample2[1][2])
                                    # Parameter errors will be written as they are 
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                    data_frame["min_BBkT_l"][i] = (get_covar_results().parmins[0])
                                    data_frame["max_BBkt_l"][i] = (get_covar_results().parmaxes[0])
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                    data_frame["min_BBNorm_l"][i] = (get_covar_results().parmins[1])
                                    data_frame["max_BBNorm_l"][i] = (get_covar_results().parmaxes[1])
                                    data_frame["BB2_kT_l"][i] = (get_covar_results().parvals[2])
                                    data_frame["min_2BBkT_l"][i] = (get_covar_results().parmins[2])
                                    data_frame["max_2BBkt_l"][i] = (get_covar_results().parmaxes[2])
                                    data_frame["BB2_Norm_l"][i] = (get_covar_results().parvals[3])
                                    data_frame["min_2BBNorm_l"][i] = (get_covar_results().parmins[3])
                                    data_frame["max_2BBNorm_l"][i] = (get_covar_results().parmaxes[3])
                            
                                    # Now the Bolometric Fluxes :
                                    data_frame["Bol_BBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0])**4.0)*get_covar_results().parvals[1])
                                    data_frame["max_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]+get_covar_results().parmaxes[0])**4.0)*(get_covar_results().parvals[1]+get_covar_results().parmaxes[1]))
                                    data_frame["min_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]-get_covar_results().parmins[0])**4.0)*(get_covar_results().parvals[1]-get_covar_results().parmins[1]))
                                    data_frame["Bol_2BBF_l"][i] = (1.076e-11*((get_covar_results().parvals[2])**4.0)*get_covar_results().parvals[3])
                                    data_frame["max_2BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[2]+get_covar_results().parmaxes[2])**4.0)*(get_covar_results().parvals[3]+get_covar_results().parmaxes[3]))
                                    data_frame["min_2BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[2]-get_covar_results().parmins[2])**4.0)*(get_covar_results().parvals[3]-get_covar_results().parmins[3]))                
                                    
                            data_frame["RStat_l"][i] = (chi)
                            data_frame["dof_l"][i] = (dof)
    
                            # For Plotting Purposes
    
                            plot_data()
                            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.05
                            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.05)
                            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                            
                            create_plots(["dbb", "sbb", "bb", "bb2"], x_max, x_min, ymax, ymin)
                            print('Finished, going to the next spectrum')
                            continue
                    if (fit_method == '4') and (get_fit_results().rstat < 1.5):
                            covar()
                            chi = get_fit_results().statval
                            dof = get_fit_results().dof
                            data_frame["RStat_l"][i] = (chi)
                            data_frame["dof_l"][i] = (dof)
    
                            parvals = np.array(get_covar_results().parvals)
                            parnames = np.array(get_covar_results().parnames)
                            parmins = np.array(get_covar_results().parmins)
                            parmaxes = np.array(get_covar_results().parmaxes)
                            #BB_kT_l[i] = (get_covar_results().parvals[0])
                            #min_BBkT_l[i] = (get_covar_results().parmins[0])
                            #max_BBkt_l[i] = (get_covar_results().parmaxes[0])
                            #BB_Norm_l[i] = (get_covar_results().parvals[1])
                            #min_BBNorm_l[i] = (get_covar_results().parmins[1])
                            #max_BBNorm_l[i] = (get_covar_results().parmaxes[1])
                    
                            # now the model unabsorbed fluxes :
                            covar()
                            cparmins = np.array(get_covar_results().parmins)
                            cparmaxes = np.array(get_covar_results().parmaxes)
                            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
    
                                    print('It seems like you have unconstrained parameters we cant calculate errors')
                                    print('No flux will be reported')
                                    
                                    # Parameter Errors will be written as 0:
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])                                
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                    
                            if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                                    print('The parameters are constrained well calculating errors')
                                    matrix = get_covar_results().extra_output
                                    is_all_zero = np. all((matrix >= 0))
                                    if is_all_zero: 
                                            sample1=sample_flux(bb,0.5,10.0, num=100, correlated=True,confidence=68)
                                            data_frame["BB_flux_l"][i] = (sample1[1][0]) 
                                            data_frame["max_BBflux_l"][i] = (sample1[1][1]-sample1[1][0])
                                            data_frame["min_BBflux_l"][i] = (sample1[1][0]-sample1[1][2])
                                            # Parameter errors will be written as they are :
                                            data_frame["BB_kT_l"][i] = (get_covar_results().parvals[1])
                                            data_frame["min_BBkT_l"][i] = (get_covar_results().parmins[1])
                                            data_frame["max_BBkt_l"][i] = (get_covar_results().parmaxes[1])
                                            data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[2])
                                            data_frame["min_BBNorm_l"][i] = (get_covar_results().parmins[2])
                                            data_frame["max_BBNorm_l"][i] = (get_covar_results().parmaxes[2])
                                            # Now the Bolometric Fluxes :
                                            data_frame["Bol_BBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0])**4.0)*get_covar_results().parvals[1])
                                            data_frame["max_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]+get_covar_results().parmaxes[0])**4.0)*(get_covar_results().parvals[1]+get_covar_results().parmaxes[1]))
                                            data_frame["min_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]-get_covar_results().parmins[0])**4.0)*(get_covar_results().parvals[1]-get_covar_results().parmins[1]))
                            plot_data()
                            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.05
                            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.05)
                            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                            create_plots(["dbb", "sbb", "bb"], x_max, x_min, ymax, ymin)
                    if (fit_method == '3'):
                            fit()
                            covar()
                            chi = get_fit_results().statval
                            dof = get_fit_results().dof
                            parvals = np.array(get_covar_results().parvals)
                            parnames = np.array(get_covar_results().parnames)
                            parmins = np.array(get_covar_results().parmins)
                            parmaxes = np.array(get_covar_results().parmaxes)
                    
                            data_frame["RStat_l"][i] = (chi)
                            data_frame["dof_l"][i] = (dof)
    
                            # now the model unabsorbed fluxes :
                            covar()
                            cparmins = np.array(get_covar_results().parmins)
                            cparmaxes = np.array(get_covar_results().parmaxes)
                            if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                                    #print('The parameters are constrained well calculating errors')
                                    #matrix = get_covar_results().extra_output
                                    #is_all_zero = np. all((matrix >= 0))
                                    #if is_all_zero: 
                                sample1=sample_flux(bb,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                                data_frame["BB_flux_l"][i] = (sample1[1][0]) 
                                data_frame["max_BBflux_l"][i] = (sample1[1][1]-sample1[1][0])
                                data_frame["min_BBflux_l"][i] = (sample1[1][0]-sample1[1][2])
                                # Parameter errors will be written as they are :
                                data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                data_frame["min_BBkT_l"][i] = (get_covar_results().parmins[0])
                                data_frame["max_BBkt_l"][i] = (get_covar_results().parmaxes[0])
                                data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                data_frame["min_BBNorm_l"][i] = (get_covar_results().parmins[1])
                                data_frame["max_BBNorm_l"][i] = (get_covar_results().parmaxes[1])
                                data_frame["fa_l"][i] = (get_covar_results().parvals[2])
                                data_frame["min_fa_l"][i] = (get_covar_results().parmins[2])
                                data_frame["max_fa_l"][i] = (get_covar_results().parmaxes[2])
                                # in this method these parameters should not be defined :
                                #data_frame["BB2_kT_l"][i] = (get_covar_results().parvals[3])
                                #data_frame["min_2BBkT_l"][i] = (get_covar_results().parmins[3])
                                #data_frame["max_2BBkt_l"][i] = (get_covar_results().parmaxes[3])
                                #data_frame["BB2_Norm_l"][i] = (get_covar_results().parvals[4])
                                #data_frame["min_2BBNorm_l"][i] = (get_covar_results().parmins[4])
                                #data_frame["max_2BBNorm_l"][i] = (get_covar_results().parmaxes[4])
                                #dbb_kT_l[i] = (get_covar_results().parvals[5])
                                #min_dbbkT_l[i] = (get_covar_results().parmins[5])
                                #max_dbbkT_l[i] = (get_covar_results().parmaxes[5])
                                #dbb_Norm_l[i] = (get_covar_results().parvals[6])
                                #min_dbbNorm_l[i] = (get_covar_results().parmins[6])
                                #max_dbbNorm_l[i] = (get_covar_results().parmaxes[6])
                            
                                # Now the Bolometric Fluxes :
                                sample_bol=sample_flux(bb,0.01,200.0, num=100, correlated=False,confidence=68)
                                data_frame["Bol_BBF_l"][i] = sample_bol[1][0]
                                data_frame["max_BolBBF_l"][i] = sample_bol[1][1]-sample_bol[1][0]
                                data_frame["min_BolBBF_l"][i] = sample_bol[1][0]-sample_bol[1][2]
                            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
    
                                    print('It seems like you have unconstrained parameters we cant calculate errors')
                                    print('No flux will be reported')
                                    
                                    # Parameter Errors will be written as 0
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                    data_frame["fa_l"][i] = (get_covar_results().parvals[2])
                                    #BB2_kT_l[i] = (get_covar_results().parvals[3])
                                    #min_2BBkT_l[i] = (0)
                                    #max_2BBkt_l[i] = (0)
                                    #BB2_Norm_l[i] = (get_covar_results().parvals[4])
                                    #min_2BBNorm_l[i] = (0)
                                    #max_2BBNorm_l[i] = (0)
                                    #dbb_kT_l[i] = (get_covar_results().parvals[5])
                                    #min_dbbkT_l[i] = (0)
                                    #max_dbbkT_l[i] = (0)
                                    #dbb_Norm_l[i] = (get_covar_results().parvals[6])
                                    #min_dbbNorm_l[i] = (0)
                                    #max_dbbNorm_l[i] = (0)
                            
                                    
                            plot_data()
                            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.1
                            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.1)
                            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                            create_plots(["fa_dbb", "fa_sbb", "bb"], x_max, x_min, ymax, ymin)
                            
                    if (fit_method == '2'):
                            covar()
                            chi = get_fit_results().statval
                            dof = get_fit_results().dof
                            parvals = np.array(get_covar_results().parvals)
                            parnames = np.array(get_covar_results().parnames)
                            parmins = np.array(get_covar_results().parmins)
                            parmaxes = np.array(get_covar_results().parmaxes)
                            #BB_kT_l[i] = (get_covar_results().parvals[0])
                            #min_BBkT_l[i] = (get_covar_results().parmins[0])
                            #max_BBkt_l[i] = (get_covar_results().parmaxes[0])
                            #BB_Norm_l[i] = (get_covar_results().parvals[1])
                            #min_BBNorm_l[i] = (get_covar_results().parmins[1])
                            #max_BBNorm_l[i] = (get_covar_results().parmaxes[1])
                            #dbb_kT_l[i] = (get_covar_results().parvals[2])
                            #min_dbbkT_l[i] = (get_covar_results().parmins[2])
                            #max_dbbkT_l[i] = (get_covar_results().parmaxes[2])
                            #dbb_Norm_l[i] = (get_covar_results().parvals[3])
                            #min_dbbNorm_l[i] = (get_covar_results().parmins[3])
                            #max_dbbNorm_l[i] = (get_covar_results().parmaxes[3])
                            #sBB_kT_l[i] = (get_covar_results().parvals[4])
                            #min_sBBkT_l[i] = (get_covar_results().parmins[4])
                            #max_sBBkt_l[i] = (get_covar_results().parmaxes[4])
                            #sBB_Norm_l[i] = (get_covar_results().parvals[5])
                            #min_sBBNorm_l[i] = (get_covar_results().parmins[5])
                            #max_sBBNorm_l[i] = (get_covar_results().parmaxes[5])
                    
                            covar()
                            cparmins = np.array(get_covar_results().parmins)
                            cparmaxes = np.array(get_covar_results().parmaxes)
                            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
                                    print('It seems like you have unconstrained parameters we cant calculate errors')
                                    print('No flux will be reported')
                                    
                                    # Parameter Errors will be written as 0:
                                    data_frame["dbb_Norm_l"][i] = (get_covar_results().parvals[1])
                                    data_frame["sBB_kT_l"][i] = (get_covar_results().parvals[2])
                                    data_frame["sBB_Norm_l"][i] = (get_covar_results().parvals[3])
                                    
    
                            if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                                    print('The parameters are constrained well calculating errors')
                                    #matrix = get_covar_results().extra_output
                                    #is_all_zero = np. all((matrix >= 0))
                                    #if is_all_zero: 
                                    sample1=sample_flux(bb,0.5,10.0, num=100, correlated=True,confidence=68)
                                    data_frame["BB_flux_l"][i] = (sample1[1][0]) 
                                    data_frame["max_BBflux_l"][i] = (sample1[1][1]-sample1[1][0])
                                    data_frame["min_BBflux_l"][i] = (sample1[1][0]-sample1[1][2])
                                    sample2=sample_flux(dbb,0.5,10.0, num=100, correlated=True,confidence=68)
                                    data_frame["dbb_flux_l"][i] = (sample2[1][0]) 
                                    data_frame["max_dbbflux_l"][i] = (sample2[1][1]-sample2[1][0])
                                    data_frame["min_dbbflux_l"][i] = (sample2[1][0]-sample2[1][2])
                                    sample3=sample_flux(sbb,0.5,10.0, num=100, correlated=True,confidence=68)
                                    data_frame["sBB_flux_l"][i] = (sample3[1][0]) 
                                    data_frame["max_sBBflux_l"][i] = (sample3[1][1]-sample3[1][0])
                                    data_frame["min_sBBflux_l"][i] = (sample3[1][0]-sample3[1][2])
                                    # Parameter errors will be written as they are :
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                    data_frame["min_BBkT_l"][i] = (get_covar_results().parmins[0])
                                    data_frame["max_BBkt_l"][i] = (get_covar_results().parmaxes[0])
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                    data_frame["min_BBNorm_l"][i] = (get_covar_results().parmins[1])
                                    data_frame["max_BBNorm_l"][i] = (get_covar_results().parmaxes[1])
                                    data_frame["dbb_kT_l"][i] = (get_covar_results().parvals[2])
                                    data_frame["min_dbbkT_l"][i] = (get_covar_results().parmins[2])
                                    data_frame["max_dbbkT_l"][i] = (get_covar_results().parmaxes[2])
                                    data_frame["dbb_Norm_l"][i] = (get_covar_results().parvals[3])
                                    data_frame["min_dbbNorm_l"][i] = (get_covar_results().parmins[3])
                                    data_frame["max_dbbNorm_l"][i] = (get_covar_results().parmaxes[3])
                                    data_frame["sBB_kT_l"][i] = (get_covar_results().parvals[4])
                                    data_frame["min_sBBkT_l"][i] = (get_covar_results().parmins[4])
                                    data_frame["max_sBBkt_l"][i] = (get_covar_results().parmaxes[4])
                                    data_frame["sBB_Norm_l"][i] = (get_covar_results().parvals[5])
                                    data_frame["min_sBBNorm_l"][i] = (get_covar_results().parmins[5])
                                    data_frame["max_sBBNorm_l"][i] = (get_covar_results().parmaxes[5])
                            
                                    #dbb_Norm_l[i] = (get_covar_results().parvals[1])
                                    #min_dbbNorm_l[i] = (get_covar_results().parmins[1])
                                    #max_dbbNorm_l[i] = (get_covar_results().parmaxes[1])
                                    #sBB_kT_l[i] = (get_covar_results().parvals[2])
                                    #min_sBBkT_l[i] = (get_covar_results().parmins[2])
                                    #max_sBBkt_l[i] = (get_covar_results().parmaxes[2])
                                    #sBB_Norm_l[i] = (get_covar_results().parvals[3])
                                    #min_sBBNorm_l[i] = (get_covar_results().parmins[3])
                                    #max_sBBNorm_l[i] = (get_covar_results().parmaxes[3])
                                    
                                    # Now the Bolometric Fluxes :
                                    data_frame["Bol_BBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0])**4.0)*get_covar_results().parvals[1])
                                    data_frame["max_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]+get_covar_results().parmaxes[0])**4.0)*(get_covar_results().parvals[1]+get_covar_results().parmaxes[1]))
                                    data_frame["min_BolBBF_l"][i] = (1.076e-11*((get_covar_results().parvals[0]-get_covar_results().parmins[0])**4.0)*(get_covar_results().parvals[1]-get_covar_results().parmins[1]))
    
                            data_frame["RStat_l"][i] = (chi)
                            data_frame["dof_l"][i] = (dof)
                            
                            
                            plot_data()
                            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.1
                            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.1)
                            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                            create_plots(["dbb", "sbb", "bb"], x_max, x_min, ymax, ymin)
                    
                    if (fit_method == '1'):
                            fit()
                            covar()
                            chi = get_fit_results().statval
                            dof = get_fit_results().dof
                            parvals = np.array(get_covar_results().parvals)
                            parnames = np.array(get_covar_results().parnames)
                            parmins = np.array(get_covar_results().parmins)
                            parmaxes = np.array(get_covar_results().parmaxes)
                               
                            #BB_kT_l[i] = (0)
                            #min_BBkT_l[i] = (0)
                            #max_BBkt_l[i] = (0)
                            #BB_Norm_l[i] = (0)
                            #min_BBNorm_l[i] = (0)
                            #max_BBNorm_l[i] = (0)
                            #BB_flux_l[i] = (0)
                            #max_BBflux_l[i] = (0)
                            #min_BBflux_l[i] = (0)                
                            
                            # now the model unabsorbed fluxes :
                            covar()
                            cparmins = np.array(get_covar_results().parmins)
                            cparmaxes = np.array(get_covar_results().parmaxes)
                            if (None in cparmins) == True or (None in cparmaxes) == True or (0 in cparmaxes) == True or (0 in cparmins) == True:
    
                                    print('It seems like you have unconstrained parameters we cant calculate errors')
                                    print('No flux will be reported')
                                    
                                    # Parameter Errors will be written as 0:
                                    data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                    data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                    #BB_flux_l[i] = (get_covar_results().parvals[2])
                                    #max_BBflux_l[i] = (0)
                                    #min_BBflux_l[i] = (0)
                            
                            else: #if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
                                print('The parameters are constrained well calculating errors')
                                #    matrix = get_covar_results().extra_output
                                #    is_all_zero = np. all((matrix >= 0))
                                #    if is_all_zero: 
                                sample1=sample_flux(bb,float(mine),float(maxe), num=100, correlated=True,confidence=68)
                                data_frame["BB_flux_l"][i] = (sample1[1][0]) 
                                data_frame["max_BBflux_l"][i] = (sample1[1][1]-sample1[1][0])
                                data_frame["min_BBflux_l"][i] = (sample1[1][0]-sample1[1][2])
                                # Parameter errors will be written as they are :
                                data_frame["BB_kT_l"][i] = (get_covar_results().parvals[0])
                                data_frame["min_BBkT_l"][i] = (get_covar_results().parmins[0])
                                data_frame["max_BBkt_l"][i] = (get_covar_results().parmaxes[0])
                                data_frame["BB_Norm_l"][i] = (get_covar_results().parvals[1])
                                data_frame["min_BBNorm_l"][i] = (get_covar_results().parmins[1])
                                data_frame["max_BBNorm_l"][i] = (get_covar_results().parmaxes[1])
                                    
                                            
                                # Now the Bolometric Fluxes :
                                    
                                sample_bol=sample_flux(bb,0.01,200.0, num=100, correlated=True,confidence=68)
                                data_frame["Bol_BBF_l"][i] = sample_bol[1][0]
                                data_frame["max_BolBBF_l"][i] = sample_bol[1][1]-sample_bol[1][0]
                                data_frame["min_BolBBF_l"][i] = sample_bol[1][0]-sample_bol[1][2]
                            
                            data_frame["RStat_l"][i] = (chi)
                            data_frame["dof_l"][i] = (dof)
                            
                            plot_data()
                            x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.1
                            x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.1)
                            ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
                            ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
                            create_plots(["dbb", "sbb", "bb"], x_max, x_min, ymax, ymin)

print_all_lengths(data_frame)

df=pd.DataFrame.from_dict(data_frame)
df.transpose()
#df.to_excel(burst_folder+source_name+'_sp_res_'+bid+'_'+fit_method+'_pow.xlsx') 
#df.to_csv(burst_folder+source_name+'_sp_res_'+bid+fit_method+'_pow.csv') 

sorted_df = df.sort_values(by='MJD_OBS_l')
sorted_df.to_excel(burst_folder+source_name+'_sp_res_'+bid+'_'+fit_method+'_pow.xlsx') 
sorted_df.to_csv(burst_folder+source_name+'_sp_res_'+bid+fit_method+'_pow.csv') 
                

    
