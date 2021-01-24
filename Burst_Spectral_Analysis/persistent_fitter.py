import glob as gl
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

# plotting the spectra, must be corrected :
def create_plots(components, x_max, x_min, ymax, ymin,index):
    print('Now Plotting the Spectrum')
    plot_fit_delchi(1,clearwindow=True, color='Black')
    fig=plt.gcf()
    ax1,ax2=fig.axes
    ax1.set_title(source_name+' BID:'+str(bid)+' MJD:'+str(mjdobs)+' Model:'+str(components))
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax2.set_xlabel('Energy [keV]', fontsize=14)
    ax1.set_ylabel('Counts/sec/keV', fontsize=14)
    ax2.set_ylabel('Sigma', fontsize=14)
    ax1.set_xlim(x_min,x_max)
    ax2.set_xlim(x_min,x_max)
    plt.savefig(burst_folder+'_b'+str(bid)+'_m'+str(index)+'_full.pdf',orientation='landscape', papertype='a4')
    plt.savefig(burst_folder+'_b'+str(bid)+'_m'+str(index)+'_full.png',orientation='landscape', papertype='a4')
    plot_fit(1,clearwindow=True,xlog=True,ylog=True, color='Black')
    if "db1" in components:
        plot_model_component("nh*db1", replot=False, overplot=True, color='Green')
    if "po1" in components:
        plot_model_component("nh*po1", replot=False, overplot=True, color='Red')
    if "bb1" in components:
        plot_model_component("nh*bb1", replot=False, overplot=True, color='Black')
    if "g1" in components:
        plot_model_component("nh*g1", replot=False, overplot=True, color='Blue')
    if "comp1" in components:        
        plot_model_component("nh*comp1", replot=False, overplot=True, color='Orange')
    plt.title(source_name+' BID:'+str(bid)+' MJD:'+str(mjdobs)+' Model:'+str(components))
    plt.xlabel('Energy [keV]', fontsize=14)
    plt.ylabel('Counts/sec/keV', fontsize=14)
    plt.xlim(x_min,x_max)
    plt.ylim(ymin,ymax)
    plt.savefig(burst_folder+'_b'+str(bid)+'_m'+str(index)+'_full_comp.pdf',orientation='landscape', papertype='a4')
    plt.savefig(burst_folder+'_b'+str(bid)+'_m'+str(index)+'_full_comp.png',orientation='landscape', papertype='a4')
    print('Spectra Plotted and files are created')

# first the constants :
source_name = 'Aql_X-1'
# enter here the burstid of the burst you would like to fit :
bid = '22'
mine=input('Enter the minimum energy of the fits :')
mines = ":"+mine
maxe=input('Enter the maximum energy of the fits :')
maxes = maxe+":"
folder = '/home/hea/ownCloud/burst_characterization_v4/'
sfolder = '/home/hea/ownCloud/burst_characterization_v4/scripts/'
burst_folder=folder+source_name+'/burst'+bid+'/'
#pers_folder = burst_folder+'/pers_analysis/'
#bkgfile = gl.glob(pers_folder+'*3c50*.pha.pi')
bkgfile = gl.glob(burst_folder+'*3c50*.pha.pi')
pre_post='pre_pers.pha'
set_stat("chi2xspecvar")
set_covar_opt("sigma",1.0)
set_conf_opt('numcores', 20)
set_conf_opt("max_rstat",250.0)
set_covar_opt('sigma',1.0)

# this script will simply fit the persistent emission of the source via several models and try list them in a table.
'''
(xstbabs.nh * (xsdiskbb.db1 + xspowerlaw.p1))
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   nh.nH        thawed            1            0       100000 10^22 atoms / cm^2
   db1.Tin      thawed            1            0         1000        keV
   db1.norm     thawed            1            0        1e+24
   po1.PhoIndex  thawed            1           -2            9
   po1.norm      thawed            1            0        1e+24
   
   (xstbabs.nh * ((xsdiskbb.db1 + xspowerlaw.p1) + xsgaussian.g1))
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   nh.nH        thawed            1            0       100000 10^22 atoms / cm^2
   db1.Tin      thawed            1            0         1000        keV
   db1.norm     thawed            1            0        1e+24
   po1.PhoIndex  thawed            1           -2            9
   po1.norm      thawed            1            0        1e+24
   g1.LineE     thawed          6.5            0        1e+06        keV
   g1.Sigma     thawed          0.1            0           10        keV
   g1.norm      thawed            1            0        1e+24

(xstbabs.nh * (xsdiskbb.db1 + xsbbodyrad.bb1))
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   nh.nH        thawed            1            0       100000 10^22 atoms / cm^2
   db1.Tin      thawed            1            0         1000        keV
   db1.norm     thawed            1            0        1e+24
   bb1.kT       thawed            3        0.001          100        keV
   bb1.norm     thawed            1            0        1e+24

(xstbabs.nh * ((xsdiskbb.db1 + xsbbodyrad.bb1) + xsgaussian.g1))
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   nh.nH        thawed            1            0       100000 10^22 atoms / cm^2
   db1.Tin      thawed            1            0         1000        keV
   db1.norm     thawed            1            0        1e+24
   bb1.kT       thawed            3        0.001          100        keV
   bb1.norm     thawed            1            0        1e+24
   g1.LineE     thawed          6.5            0        1e+06        keV
   g1.Sigma     thawed          0.1            0           10        keV
   g1.norm      thawed            1            0        1e+24

(xstbabs.nh * xscomptt.comp1)
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   nh.nH        thawed            1            0       100000 10^22 atoms / cm^2
   comp1.redshift frozen            0       -0.999           10
   comp1.T0     thawed          0.1         0.01          100        keV
   comp1.kT     thawed           50            2          500        keV
   comp1.taup   thawed            1         0.01          100
   comp1.approx frozen            1            0            5
   comp1.norm   thawed            1            0        1e+24
   '''
# the models to be fit, these are also individual cases : 

models = ['xstbabs.nh*(xsdiskbb.db1+xspowerlaw.po1)',
          'xstbabs.nh*(xsdiskbb.db1+xspowerlaw.po1+xsgaussian.g1)',
          'xstbabs.nh*(xsdiskbb.db1+xsbbodyrad.bb1)','xstbabs.nh*(xsdiskbb.db1+xsbbodyrad.bb1+xsgaussian.g1)'
          ,'xstbabs.nh*(xscomptt.comp1)']

# the columns of the resulting table can be : 

# Source Name, BurstID, MJD, DATEOBS, EXP, PRE/POST, model_name, NH, minNH, maxNH, Tin/T0, minTin, maxTin, dbb_norm, mindbb_norm, maxdbb_norm, 
# phoind/bb1kT/compkT, minphoind, maxphoind, pow_norm, minpow_norm, maxpow_norm, g1lineE,ming1lineE, maxg1lineE, g1sigma, ming1sigma, maxg1sigma, g1norm, ming1norm, maxg1norm, 
# comptau, mincomptau, maxcomptau, compnorm, mincompnorm, maxcompnorm, total_flux, mintotal_flux, max_total_flux, mine, maxe, chi2, dof, rchi

col_names = ["Source Name", "BurstID", "MJD", "DATEOBS", "EXP",\
             "PRE/POST", "model_name", "NH", "minNH", "maxNH","Tin/T0",\
             "minTin", "maxTin", "dbb_norm", "mindbb_norm", \
             "maxdbb_norm", "phoind/bb1kT/compkT", "minphoind", "maxphoind", \
             "pow_norm", "minpow_norm", "maxpow_norm", "g1lineE", "ming1lineE",\
             "maxg1lineE", "g1sigma", "ming1sigma", "maxg1sigma",\
             "g1norm", "ming1norm", "maxg1norm", "comptau",\
             "mincomptau", "maxcomptau", "total_flux",\
             "mintotal_flux", "maxtotal_flux", "mine", "maxe","chi2", "dof", "rchi2"]

# definition of the dictionary to keep the best fit parameters :

data_frame = {col: [0]*len(models) for col in col_names}

for i in range(len(models)):
    # first lets try the pre spectrum : 
    spec = gl.glob(burst_folder+pre_post)
    sp_hdu = fits.open(spec[0])
    print('read src spec: '+str(spec[0]))
    mjdobs = sp_hdu[1].header['MJD-OBS']
    date_obsi = sp_hdu[1].header['DATE-OBS']
    exposure = sp_hdu[1].header['EXPOSURE']
    obsid = sp_hdu[1].header['OBS_ID']
    date_obs = str(Time(date_obsi,format='isot', scale='utc'))
    data_frame["Source Name"][i]=source_name
    data_frame["BurstID"][i] = bid
    data_frame["MJD"][i] = mjdobs
    data_frame["DATEOBS"][i] = date_obs
    data_frame["EXP"][i]=exposure
    data_frame["mine"][i] = mine
    data_frame["maxe"][i] = maxe
    data_frame["PRE/POST"][i] = pre_post
    data_frame["model_name"][i]=models[i]
        
    # first for the pre burst spectrum :
    load_pha(1, str(spec[0]),use_errors=True)
    load_arf(1, sfolder+'nixtionaxis20170601_combined_v004_1434.arf')
    load_rmf(1, sfolder+'nixti20170601_combined_v002_1434.rmf')
    #load_arf(1, sfolder+'nixtionaxis20170601_combined_v004.arf')
    #load_rmf(1, sfolder+'nixti20170601_combined_v002.rmf')
    load_bkg(1, bkgfile[0],use_errors=True)
    subtract()
    print('This script only subtracts ni3C50 background')
    print('Ignoring : '+mines+' '+maxes)
    ignore(mines+','+maxes)
    print('Grouping the data to have at least 50 counts per channel')
    group_counts(1, 50)
    plot_data()
    x_max=max(get_data_plot().x)+max(get_data_plot().x)*0.05
    x_min = np.abs(min(get_data_plot().x)-min(get_data_plot().x)*0.05)
    ymax = max(get_data_plot().y)+max(get_data_plot().y)*0.2
    ymin = np.abs(min(get_data_plot().y)-min(get_data_plot().y)*0.05)
    set_source(models[i])
    if i == 4:
        comp1.redshift=0.0
        freeze(comp1.redshift)
        comp1.approx=1.1
        freeze(comp1.approx)
        thaw(comp1.T0)
        thaw(comp1.kT)
        thaw(comp1.norm)
        comp1.norm=1e-2
    if ((i == 1) or (i == 3)):
        g1.LineE = 6.4
        g1.Sigma.min=0.01
        g1.Sigma.max=0.7
        g1.norm=1e-3
        g1.LineE.min=6.0
        g1.LineE.max=7.0
    nh.nH=1.0
    freeze(nh.nh)
    db1.Tin=0.8
    freeze(db1.Tin)
    fit()
    thaw(nh.nH)
    thaw(db1.Tin)
    set_method('moncar')
    fit()
    set_method('levmar')
    fit()
    conf()
    # we must record the outputs here withe something like the following :
    chi = get_fit_results().statval
    dof = get_fit_results().dof
    parvals = np.array(get_conf_results().parvals)
    parnames = np.array(get_conf_results().parnames)
    parmins = np.array(get_conf_results().parmins)
    parmaxes = np.array(get_conf_results().parmaxes)
    data_frame['chi2'][i] = chi
    data_frame['dof'][i] = dof
    data_frame['rchi2'][i] = chi/dof
    data_frame['NH'][i] = parvals[0]
    data_frame['minNH'][i] = parmins[0]
    data_frame['maxNH'][i] = parmaxes[0]
    data_frame["Tin/T0"][i] = parvals[1]
    data_frame["minTin"][i] = parmins[1]
    data_frame["maxTin"][i] = parmaxes[1]
    data_frame["dbb_norm"][i]= parvals[2]
    data_frame["mindbb_norm"][i] = parmins[2]
    data_frame["maxdbb_norm"][i] = parmaxes[2]
    data_frame["phoind/bb1kT/compkT"][i] = parvals[3]
    data_frame["minphoind"][i] = parmins[3]
    data_frame["maxphoind"][i] = parmaxes[3]
    data_frame["pow_norm"][i] = parvals[4]
    data_frame["minpow_norm"][i] = parmins[4]
    data_frame['maxpow_norm'][i] = parmaxes[4]
    if ((i == 1) or (i == 3)):
        data_frame['g1lineE'][i] = parvals[5]
        data_frame['ming1lineE'][i] = parmins[5]
        data_frame['maxg1lineE'][i] = parmaxes[5]
        data_frame['g1sigma'][i] = parvals[6]
        data_frame['ming1sigma'][i] = parmins[6]
        data_frame['maxg1sigma'][i] = parmaxes[6]
        data_frame['g1norm'][i]=parvals[7]
        data_frame['ming1norm'][i]=parmins[7]
        data_frame['maxg1norm'][i]=parmaxes[7]    
    if i == 4:    
        data_frame["Tin/T0"][i] = parvals[1]
        data_frame["minTin"][i] = parmins[1]
        data_frame["maxTin"][i] = parmaxes[1]
        data_frame["phoind/bb1kT/compkT"][i] = parvals[2]
        data_frame["minphoind"][i] = parmins[2]
        data_frame["maxphoind"][i] = parmaxes[2]
        data_frame["comptau"][i] = parvals[3]
        data_frame["mincomptau"][i] = parmins[3]
        data_frame["maxcomptau"][i] = parmaxes[3]
        data_frame["pow_norm"][i] = parvals[4]
        data_frame["minpow_norm"][i] = parmins[4]
        data_frame["maxpow_norm"][i] = parmaxes[4]
    # now the flux :
    if i == 4:
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)
        if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
            try:
                sample_bol=sample_flux(comp1,0.01,200.0, num=1000, correlated=True,confidence=68)
                data_frame["total_flux"][i] = sample_bol[1][0]
                data_frame["mintotal_flux"][i] = sample_bol[1][1]-sample_bol[1][0]
                data_frame["maxtotal_flux"][i] = sample_bol[1][0]-sample_bol[1][2]
            except:
                data_frame["total_flux"][i] = 0
                data_frame["mintotal_flux"][i] = 0
                data_frame["maxtotal_flux"][i] = 0
        create_plots(["comp1"], x_max, x_min, ymax, ymin,i)
    if i == 0:
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)
        if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
            try:
                sample_bol=sample_flux(db1+po1,0.01,200.0, num=1000, correlated=True,confidence=68)
                data_frame["total_flux"][i] = sample_bol[1][0]
                data_frame["mintotal_flux"][i] = sample_bol[1][1]-sample_bol[1][0]
                data_frame["maxtotal_flux"][i] = sample_bol[1][0]-sample_bol[1][2]
            except:
                data_frame["total_flux"][i] = 0
                data_frame["mintotal_flux"][i] = 0
                data_frame["maxtotal_flux"][i] = 0
        create_plots(["dbb","po1"],x_max,x_min,ymax, ymin,i)
    if i == 1:
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)
        if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
            try:
                sample_bol=sample_flux(db1+po1+g1,0.01,200.0, num=1000, correlated=True,confidence=68)
                data_frame["total_flux"][i] = sample_bol[1][0]
                data_frame["mintotal_flux"][i] = sample_bol[1][1]-sample_bol[1][0]
                data_frame["maxtotal_flux"][i] = sample_bol[1][0]-sample_bol[1][2]
            except:
                data_frame["total_flux"][i] = 0
                data_frame["mintotal_flux"][i] = 0
                data_frame["maxtotal_flux"][i] = 0
        create_plots(["dbb","po1","g1"],x_max,x_min,ymax, ymin,i)
    if i == 2:
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)
        if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
            try:
                sample_bol=sample_flux(db1+bb1,0.01,200.0, num=1000, correlated=True,confidence=68)
                data_frame["total_flux"][i] = sample_bol[1][0]
                data_frame["mintotal_flux"][i] = sample_bol[1][1]-sample_bol[1][0]
                data_frame["maxtotal_flux"][i] = sample_bol[1][0]-sample_bol[1][2]
            except:
                data_frame["total_flux"][i] = 0
                data_frame["mintotal_flux"][i] = 0
                data_frame["maxtotal_flux"][i] = 0
        create_plots(["dbb","bb1"],x_max,x_min,ymax, ymin,i)
    if i == 3:
        covar()
        cparmins = np.array(get_covar_results().parmins)
        cparmaxes = np.array(get_covar_results().parmaxes)
        if (None in cparmins) == False and (None in cparmaxes) == False and (0 in cparmaxes) == False and (0 in cparmins) == False:
            try:
                sample_bol=sample_flux(db1+bb1+g1,0.01,200.0, num=1000, correlated=True,confidence=68)
                data_frame["total_flux"][i] = sample_bol[1][0]
                data_frame["mintotal_flux"][i] = sample_bol[1][1]-sample_bol[1][0]
                data_frame["maxtotal_flux"][i] = sample_bol[1][0]-sample_bol[1][2]
            except:
                data_frame["total_flux"][i] = 0
                data_frame["mintotal_flux"][i] = 0
                data_frame["maxtotal_flux"][i] = 0
        create_plots(["dbb","bb1","g1"],x_max,x_min,ymax, ymin,str(i))

df=pd.DataFrame.from_dict(data_frame)
df.transpose()
df.to_excel(burst_folder+source_name+'_sp_res_'+str(bid)+'_pre_pers.xlsx') 
df.to_csv(burst_folder+source_name+'_sp_res_'+str(bid)+'_pre_pres.csv') 

        
        


