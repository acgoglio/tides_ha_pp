# Emanuel Clementi 22-02-2019
# This script plots the timeseries diagnostics for 2 experiments


import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from netCDF4 import Dataset
import netCDF4 as ncdf
import datetime
import pandas as pd
import glob
from numpy import *
import warnings
warnings.filterwarnings("ignore")

ref1_name='MFSe1r1'
ref2_name='MFSe1r1_2017'
exp1_name='EAS1r1'
exp2_name='EAS1r1_2'

pathdir_ref1='/work/md04916/exp/MFS_EAS1/cf_exp/exp_mrsp.sh/diag_base_eas1_pp.xml/' 
pathdir_ref2='/work/md04916/exp/MFS_EAS1/cf_expm/exp_mrsp.sh/diag_base_eas1_cfm.xml/'
pathdir_exp1='/work/md04916/exp/MFS_EAS1/myp_ass/exp_mrsp.sh/diag_base_eas1.xml/'
pathdir_exp2='/work/md04916/exp/MFS_EAS1/myp3_ass/exp_mrsp.sh/diag_base_eas1.xml/'
#/work/md04916/cosmo_validazione/diag_test_eas1/']

figdir='/work/ec04916/FIGURE/DIAG/Diff_'+ref1_name+'_'+exp2_name+'/'

T_avg=np.matrix([[2016,1,1],[2017,12,31]])

filename=['tra_t_gb_ts','tra_p_gb_ts','tra_n_gb_ts','tra_t_sc_ts','tra_p_sc_ts','tra_n_sc_ts','tra_t_ot_ts', \
'tra_p_ot_ts','tra_n_ot_ts','tra_t_co_ts','tra_p_co_ts','tra_n_co_ts','tra_t_me_ts','tra_p_me_ts','tra_n_me_ts',\
'T.nc_SST_ts','T.nc_0_150_ts','T.nc_150_600_ts','T.nc_600_btm_ts','T.nc_basin_ts','S.nc_SSS_ts','S.nc_0_150_ts','S.nc_150_600_ts',\
'S.nc_600_btm_ts','S.nc_basin_ts','sossheig_ts',\
#,'somxl010_ts','sohefldo_ts','sowaflup_ts','soevapor_ts','soprecip_ts',\
#'sorunoff_ts','soshfldo_ts','solofldo_ts','sosefldo_ts','solafldo_ts','velmodw_ts',
'V.nc_SSV_ts','V.nc_0_150_ts','V.nc_150_600_ts','V.nc_600_btm_ts','V.nc_basin_ts','K.nc_SSK_ts','K.nc_0_150_ts','K.nc_150_600_ts','K.nc_600_btm_ts','K.nc_basin_ts']

varname=['transptx','transppx','transpnx','transpty','transppy','transpny','transpty','transppy','transpny','transpty',\
'transppy','transpny','transpty','transppy','transpny','votemper','votemper','votemper','votemper','votemper','vosaline',\
'vosaline','vosaline','vosaline','vosaline','sossheig',\
#'somxl010','sohefldo','sowaflup','soevapor','soprecip','sorunoff',\
#'soshfldo','solofldo','sosefldo','solafldo','velmodw',i
'velmod','velmod','velmod','velmod','velmod','vokenerg','vokenerg','vokenerg','vokenerg','vokenerg']

varlongname=['Net Transport Gibraltar','Eastward Transport Gibraltar','Westward Transport Gibraltar',\
'Net Transport Sicily Channel','Northward Transport Sicily Channel','Southward Transport Sicily Channel',\
'Net Transport Otranto Channel','Northward Transport Otranto Channel','Southward Transport Otranto Channel',\
'Net Transport Corsica Channel','Northward Transport Corsica Channel','Southward Transport Corsica Channel',\
'Net Transport Messina Strait','Northward Transport Messina Strait','Southward Transport Messina Strait',\
'Sea Surface Temperature','Temperature 0-150m','Temperature 150-600m','Temperature 600m-bottom','Basin averaged Temperature',\
'Sea Surface Salinity','Salinity 0-150m','Salinity 150-600m','Salinity 600m-bottom','Basin averaged Salinity',\
'Sea Surface Height',\
#'Mixed Layer Depth','Net Downward Heat Flux','Net Upward Water Flux','Water Evaporation Flux',\
#'Precipitation Flux','River runoffs','Downward Shortwave Radiation','Downward Longwave Radiation',\
#'Downward Sensible Heat Flux','Downward Latent Heat Flux','Wind Stress',\
'Sea Surface Current','Current 0-150m','Current 150-600m','Current 600m-bottom','Basin averaged Current',\
'Sea Surface Kinetic Energy','Kinetic Energy 0-150m','Kinetic Energy 150-600m','Kinetic Energy 600m-bottom','Basin averaged Kinetic Energy']

units=['[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[$^\circ$C]',\
'[$^\circ$C]','[$^\circ$C]','[$^\circ$C]','[$^\circ$C]','[PSU]','[PSU]','[PSU]','[PSU]','[PSU]','[cm]',\
#'[m]','[W/m$^2$]',\
#'[Kg/m$^2$/s]','[Kg/m$^2$/s]','[Kg/m$^2$/s]','[Kg/m$^2$/s]',\
#'[W/m$^2$]','[W/m$^2$]','[W/m$^2$]','[W/m$^2$]','[m/s]',
'[m/s]','[m/s]','[m/s]','[m/s]','[m/s]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]']

for n in range(0, len(filename)):
# Experiment N. 1
        file_ref1=pathdir_ref1+filename[n]+'.nc'
        print(file_ref1)
        fh_ref1=ncdf.Dataset(file_ref1,'r')
        var_ref1 = fh_ref1.variables[varname[n]][:]
        aa=shape(var_ref1)
        if len(aa)==3:
          VAR_ref1=np.squeeze(var_ref1[:,0,0])
        if len(aa)==4:
          VAR_ref1=np.squeeze(var_ref1[:,0,0,0])
        time_ref1 = fh_ref1.variables['time'][:]
        fh_ref1.close()
        ndays1=len(time_ref1)
        data_ini1=datetime.date(1969, 12, 31) + datetime.timedelta(days=time_ref1[0]/(24*60*60))
        data_fin1=datetime.date(1969, 12, 31) + datetime.timedelta(days=time_ref1[-1]/(24*60*60)+1)
        print(data_ini1)
        print(data_fin1)
        VAR_ref1=np.concatenate((VAR_ref1,NaN,NaN),axis=None)
        VAR_ref1=VAR_ref1[1:-1]


# Experiment N. 2
        file_ref2=pathdir_ref2+filename[n]+'.nc'
        print(file_ref2)
        fh_ref2=ncdf.Dataset(file_ref2,'r')
        var_ref2 = fh_ref2.variables[varname[n]][:]
        aa=shape(var_ref2)
        if len(aa)==3:
          VAR_ref2=np.squeeze(var_ref2[:,0,0])
        if len(aa)==4:
          VAR_ref2=np.squeeze(var_ref2[:,0,0,0])
        time_ref2 = fh_ref2.variables['time'][:]
        fh_ref2.close()
        ndays2=len(time_ref2)
        data_ini2=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_ref2[0]/(24*60))
        data_fin2=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_ref2[-1]/(24*60)+1)
        print(data_ini2)
        print(data_fin2)

    #    idx_avg_ini=(datetime.date(T_avg[0,0],T_avg[0,1],T_avg[0,2]) - data_ini).days
    #    idx_avg_fin=(datetime.date(T_avg[1,0],T_avg[1,1],T_avg[1,2]) - data_ini).days
    #    print(idx_avg_ini)
    #    print(idx_avg_fin)
#        times_ref1 = pd.date_range(str(data_ini1),periods=ndays, freq = "1d")
        VAR_ref=np.concatenate((VAR_ref1,VAR_ref2), axis=None)
        time_ref=np.concatenate((time_ref1,time_ref2), axis=None)
        ndays=len(time_ref)
        times_ref = pd.date_range(str(data_ini1),periods=ndays, freq = "1d")
        #data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[0]/(24*60))
        #data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[-1]/(24*60)+1)
        #print(data_ini)
        #print(data_fin)

        file_exp1=pathdir_exp1+filename[n]+'.nc'
        print(file_exp1)
        fh_exp1=ncdf.Dataset(file_exp1,'r')
        var_exp1 = fh_exp1.variables[varname[n]][:]
        aa=shape(var_exp1)
        if len(aa)==3:
          VAR_exp1=np.squeeze(var_exp1[:,0,0])
        if len(aa)==4:
          VAR_exp1=np.squeeze(var_exp1[:,0,0,0])
        time_exp1 = fh_exp1.variables['time'][:]
        fh_exp1.close()
        ndays=len(time_exp1)
        data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_exp1[0]/(24*60))
        data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_exp1[-1]/(24*60)+1)
        print(data_ini)
        print(data_fin)
        times_exp1 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")

        file_exp2=pathdir_exp2+filename[n]+'.nc'
        print(file_exp2)
        fh_exp2=ncdf.Dataset(file_exp2,'r')
        var_exp2 = fh_exp2.variables[varname[n]][:]
        aa=shape(var_exp2)
        if len(aa)==3:
          VAR_exp2=np.squeeze(var_exp2[:,0,0])
        if len(aa)==4:
          VAR_exp2=np.squeeze(var_exp2[:,0,0,0])
        time_exp2 = fh_exp2.variables['time'][:]
        fh_exp2.close()
        ndays=len(time_exp2)
        data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_exp2[0]/(24*60))
        data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time_exp2[-1]/(24*60)+1)
        print(data_ini)
        print(data_fin)
        times_exp2 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")


        fig = plt.figure(n,figsize=(9,7.5))
        ax = fig.add_subplot(111)
#       ax.plot(times1,VAR1,colexp[exp],label=expname[exp]+': AVG= '+str(MEAN1)+units[n],linewidth=1.5)
        print(len(times_ref))
        print(len(VAR_ref))
        print(len(VAR_exp1))
        ax.plot(times_ref,VAR_ref-VAR_exp1,'b',label='EAS1r1',linewidth=1.5)
        ax.plot(times_ref,VAR_ref-VAR_exp2,'g',label='EAS1r1_2',linewidth=1.5)
        ax.grid('on')
        leg = plt.legend(loc='upper left', ncol=2,  shadow=True, fancybox=True, fontsize=12)
        leg.get_frame().set_alpha(0.3)
        plt.title('DIFF('+ref1_name+' - '+exp1_name+'):  '+ varlongname[n]+' '+units[n] ,fontsize=18)
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_minor_locator(mdates.MonthLocator((1,4,7,10)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("\n%Y"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
        plt.savefig(figdir+filename[n]+'.jpg')
