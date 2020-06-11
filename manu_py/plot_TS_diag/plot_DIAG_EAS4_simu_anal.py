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

expname=['simu_EAS4']
#expname=['simu_EAS4','EAS4']
colexp=['r']
#colexp=['r','b']
nexp=len(expname)

pathdir=['/work/ec04916/exp/eas5/simu_EAS4_a/exp_mrsp.sh/diag_base.xml/']
#pathdir=['/work/ec04916/exp/eas5/simu_EAS4_a/exp_mrsp.sh/diag_base.xml',\
#'/work/ec04916/exp/eas5/EAS4_a_PJ/exp_mrsp.sh/diag_base.xml']

#figdir='/work/ec04916/FIGURE/DIAG/'+expname[0]+'_'+expname[1]+'/'
figdir='/work/ec04916/FIGURE/DIAG/'+expname[0]+'/'

T_avg=np.matrix([[2015,1,1],[2015,12,31]])

filename=['tra_t_gb_ts','tra_p_gb_ts','tra_n_gb_ts','tra_t_sc_ts','tra_p_sc_ts','tra_n_sc_ts','tra_t_ot_ts', \
'tra_p_ot_ts','tra_n_ot_ts','tra_t_co_ts','tra_p_co_ts','tra_n_co_ts','tra_t_me_ts','tra_p_me_ts','tra_n_me_ts',\
'T.nc_SST_ts','T.nc_0_150_ts','T.nc_150_600_ts','T.nc_600_btm_ts','T.nc_basin_ts','S.nc_SSS_ts','S.nc_0_150_ts','S.nc_150_600_ts',\
'S.nc_600_btm_ts','S.nc_basin_ts','sossheig_ts','somxl010_ts','sohefldo_ts','sowaflup_ts','soevapor_ts','soprecip_ts',\
'sorunoff_ts','soshfldo_ts','solofldo_ts','sosefldo_ts','solafldo_ts','velmodw_ts','V.nc_SSV_ts','V.nc_0_150_ts',\
'V.nc_150_600_ts','V.nc_600_btm_ts','V.nc_basin_ts','K.nc_SSK_ts','K.nc_0_150_ts','K.nc_150_600_ts','K.nc_600_btm_ts','K.nc_basin_ts']

varname=['transptx','transppx','transpnx','transpty','transppy','transpny','transpty','transppy','transpny','transpty',\
'transppy','transpny','transpty','transppy','transpny','votemper','votemper','votemper','votemper','votemper','vosaline',\
'vosaline','vosaline','vosaline','vosaline','sossheig','somxl010','sohefldo','sowaflup','soevapor','soprecip','sorunoff',\
'soshfldo','solofldo','sosefldo','solafldo','velmodw','velmod','velmod','velmod','velmod','velmod',\
'vokenerg','vokenerg','vokenerg','vokenerg','vokenerg']

varlongname=['Net Transport Gibraltar','Eastward Transport Gibraltar','Westward Transport Gibraltar',\
'Net Transport Sicily Channel','Northward Transport Sicily Channel','Southward Transport Sicily Channel',\
'Net Transport Otranto Channel','Northward Transport Otranto Channel','Southward Transport Otranto Channel',\
'Net Transport Corsica Channel','Northward Transport Corsica Channel','Southward Transport Corsica Channel',\
'Net Transport Messina Strait','Northward Transport Messina Strait','Southward Transport Messina Strait',\
'Sea Surface Temperature','Temperature 0-150m','Temperature 150-600m','Temperature 600m-bottom','Basin averaged Temperature',\
'Sea Surface Salinity','Salinity 0-150m','Salinity 150-600m','Salinity 600m-bottom','Basin averaged Salinity',\
'Sea Surface Height','Mixed Layer Depth','Net Downward Heat Flux','Net Upward Water Flux','Water Evaporation Flux',\
'Precipitation Flux','River runoffs','Downward Shortwave Radiation','Downward Longwave Radiation',\
'Downward Sensible Heat Flux','Downward Latent Heat Flux','Wind Stress',\
'Sea Surface Current','Current 0-150m','Current 150-600m','Current 600m-bottom','Basin averaged Current',\
'Sea Surface Kinetic Energy','Kinetic Energy 0-150m','Kinetic Energy 150-600m','Kinetic Energy 600m-bottom','Basin averaged Kinetic Energy']

units=['[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[Sv]','[$^\circ$C]',\
'[$^\circ$C]','[$^\circ$C]','[$^\circ$C]','[$^\circ$C]','[PSU]','[PSU]','[PSU]','[PSU]','[PSU]','[cm]','[m]','[W/m$^2$]',\
'[Kg/m$^2$/s]','[Kg/m$^2$/s]','[Kg/m$^2$/s]','[Kg/m$^2$/s]',\
'[W/m$^2$]','[W/m$^2$]','[W/m$^2$]','[W/m$^2$]','[m/s]','[m/s]','[m/s]','[m/s]','[m/s]','[m/s]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]','[m$^2$/s$^2$]']

for exp in range(0,nexp):
   for n in range(15, len(filename)):
# Experiment N. 1
        file1=pathdir[exp]+filename[n]+'.nc'
        print(file1)
        fh1=ncdf.Dataset(file1,'r')
        var1 = fh1.variables[varname[n]][:]
        aa=shape(var1)
        if len(aa)==3:
          VAR1=np.squeeze(var1[:,0,0])
        if len(aa)==4:
          VAR1=np.squeeze(var1[:,0,0,0])
        time1 = fh1.variables['time'][:]
        fh1.close()
#Evaluates the time period
        ndays=len(time1)
        data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[0]/(24*60))
        print(data_ini)
        data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[-1]/(24*60)+1)
        print(data_fin)
        times1 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")
# Evaluates the average value inthe specified period (T_avg)        
        idx_avg_ini=(datetime.date(T_avg[0,0],T_avg[0,1],T_avg[0,2]) - data_ini).days
        idx_avg_fin=(datetime.date(T_avg[1,0],T_avg[1,1],T_avg[1,2]) - data_ini).days
        MEAN1=np.mean(VAR1[idx_avg_ini:idx_avg_fin])

        fig = plt.figure(n,figsize=(9,7.5))
        ax = fig.add_subplot(111)
        ax.plot(times1,VAR1,colexp[exp],label=expname[exp]+': AVG= '+str(MEAN1)+units[n],linewidth=1.5)
        leg = plt.legend(loc='upper left', ncol=2,  shadow=True, fancybox=True, fontsize=12)
        leg.get_frame().set_alpha(0.3)
        ax.grid('on')
        plt.title( varlongname[n]+' '+units[n] ,fontsize=18)
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_minor_locator(mdates.MonthLocator((1,4,7,10)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("\n%Y"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
        plt.savefig(figdir+filename[n]+'.jpg')

