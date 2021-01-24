# this code will read the results of plot the results of time resolved spectral analysis of X-ray bursts :

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
from scipy import integrate
import math
# Enter here the name of the source as in the burster_v3f.dat : 
source_name = '4U_1608-522'
# enter here the burstid of the burst you would like to fit :
bid = '5'
# Enter here the fit_method for which the results will be plotted :

fit_method = '1'

# maximum time range for the plot
max_time = 35.0
mine=input('Enter the minimum energy of the fits :')
mines = ":"+mine
maxe=input('Enter the maximum energy of the fits :')
maxes = ":"+maxe

burst_folder = '/Users/tolga/ownCloud/burst_characterization_v4/'
#burst_folder ='/home/hea/ownCloud/burst_characterization_v4/'
sfolder = '/Users/tolga/ownCloud/burst_characterization_v4/scripts/Burst_Spectral_Analysis/'
#sfolder = '/home/hea/ownCloud/burst_characterization_v3/scripts'

data = pd.read_csv(burst_folder+source_name+'/burst'+bid+'/'+source_name+'_sp_res_'+bid+fit_method+'_pow.csv')

Source_Name_l=data['Source_Name_l']
BID_l=data['BID_l'].to_numpy()
SID_l=data['SID_l'].to_numpy()
OBS_ID_l=data['OBS_ID_l'].to_numpy()
MJD_OBS_l=data['MJD_OBS_l'].to_numpy()
exp_l=data['exp_l'].to_numpy()
#b_time = MJD_OBS_l-MJD_OBS_l[0]*86400+exp_l/2.0
b_time=(MJD_OBS_l-MJD_OBS_l.min())*86400+exp_l/2.0-0.5
DATE_OBS_l=data['DATE_OBS_l'].to_numpy()
exp_err = [exp_l/2.0,exp_l/2.0]
NH_l=data['NH_l'].to_numpy()
BB_kT_l=data['BB_kT_l'].to_numpy()
BBkT_err=[np.abs(data['min_BBkT_l']),data['max_BBkt_l']]
BB_Norm_l=data['BB_Norm_l'].to_numpy()
min_BBNorm_l=data['min_BBNorm_l'].to_numpy()
max_BBNorm_l=data['max_BBNorm_l'].to_numpy()
min_BBkT_l=data['min_BBkT_l'].to_numpy()
max_BBkt_l=data['max_BBkt_l'].to_numpy()
BBNorm_err=[np.abs(data['min_BBNorm_l']),data['max_BBNorm_l']]
BB_flux_l=data['BB_flux_l'].to_numpy()
min_BBflux_l=data['min_BBflux_l'].to_numpy()
max_BBflux_l=data['max_BBflux_l'].to_numpy()
BB_flux_err = [data['min_BBflux_l']/1e-9,data['max_BBflux_l']/1e-9]
Bol_BBF_l=data['Bol_BBF_l'].to_numpy()/1e-9
min_BolBBF_l=data['min_BolBBF_l'].to_numpy()
max_BolBBF_l=data['max_BolBBF_l'].to_numpy()

sel_burst_time = np.where(Bol_BBF_l > max(Bol_BBF_l)*0.05)

fluence_flux = Bol_BBF_l[sel_burst_time[0]]
fluence_flux_err = ((min_BolBBF_l[sel_burst_time[0]]+max_BolBBF_l[sel_burst_time[0]])/2.0)/1e-9
fl_exp = exp_l[sel_burst_time[0]]
fluence_int = integrate.trapz(fluence_flux, x=None, dx=fl_exp[0])

sum_fl = 0
for i in range(len(fluence_flux)):
    sum_fl = sum_fl+fl_exp[i]**2.0*(fluence_flux_err[i]**2.0)

error_of_sum = math.sqrt(sum_fl)

print("Fluence :")
print((fluence_int*1e-9))
print("Error of fluence :")
print((error_of_sum*1e-9))


 