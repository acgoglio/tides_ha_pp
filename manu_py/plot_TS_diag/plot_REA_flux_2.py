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

expname=['MFSe1r1','EAS1r1_2','EAS1r1_2']
colexp=['r','g','g']
nexp=len(expname)


pathdir=['/work/md04916/exp/MFS_EAS1/cf_expm/exp_mrsp.sh/diag_base_eas1_cfm.xml/',\
'/work/md04916/exp/MFS_EAS1/myp3_ass/exp_mrsp.sh/diag_base_eas1.xml/',\
'/work/md04916/exp/MFS_EAS1/myp3_ass8/exp_mrsp.sh/diag_base_eas1.xml/']

#cf_expm = REA16 Claudia da model output
#myp3_ass = REA16 Max reanalisi_2 2016-2017 (cambio SSH boundary e MDT)
#myp3_ass8 = REA16 Max reanalisi_2 2018 (continuazione di myp3_ass) 

figdir='/work/ec04916/FIGURE/DIAG/REA_FLUXES/'


filename=['sossheig_ts','soshfldo_ts','sohefldo_ts','sowaflup_ts']

varname=['sossheig','soshfldo','sohefldo','sowaflup']

varlongname = ['Sea Surface Height','Shortwave Radiation','Net Downward Heat Flux','Net Upward Water Flux']

print(VAR.shape)
for exp in range(0,nexp):
   print(exp)
   for n in range(0, len(filename)):
        print(n)
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
        print ndays

        data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[0]/(24*60))
        data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time1[-1]/(24*60)+1)
        print(data_ini)
        print(data_fin)

        times1 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")
        #print(times1)
# Evaluates the average value inthe specified period (T_avg)        

        fig = plt.figure(n,figsize=(9,7.5))
        ax = fig.add_subplot(111)
#        ax.plot(times1,VAR1,colexp[exp],label=expname[exp]+': AVG= '+str(MEAN1)+units[n],linewidth=1.5)
#        ax.plot(times1,VAR1,colexp[exp],label=expname[exp],linewidth=1.5)
        if exp==0 or exp==1:
          ax.plot(times1,VAR1,colexp[exp],label=expname[exp],linewidth=1.5)
          leg = plt.legend(loc='upper left', ncol=2,  shadow=True, fancybox=True, fontsize=12)
          leg.get_frame().set_alpha(0.3)
        else:
          ax.plot(times1,VAR1,colexp[exp],linewidth=1.5)

        ax.grid('on')
        plt.title( varlongname[n] ,fontsize=18)
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_minor_locator(mdates.MonthLocator((1,4,7,10)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("\n%Y"))
        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
        plt.savefig(figdir+filename[n]+'.jpg')

