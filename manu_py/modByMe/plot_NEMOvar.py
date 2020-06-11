# Written by Emanuela Clementi 22-02-2019
# Modified by Anna Chiara Goglio 06-06-2019
# This script plots NEMO var timeseries 

# Setting the env
print("Setting the environment...")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import os
from netCDF4 import Dataset
import netCDF4 as ncdf
import datetime
import pandas as pd
import glob
from numpy import *
import warnings
warnings.filterwarnings("ignore")
print("..Done!")

# Exp name
expname=['EAS5']
print("Exp names: "+expname[0] )

# Input path
pathdir='/work/work/ag15419/exp/eas5/simu_ctrl0/output/201601/'
print("Inputh path: "+pathdir)

# Output path
figdir='/work/ec04916/FIGURE/DIAG/'+expname[0]+'/'
print("Output path: "+figdir)


# Time range
T_avg=np.array([[2016, 1, 1], [2016, 1, 3]])
print("Time interval: "+str(T_avg[0][2])+"/"+str(T_avg[0][1])+"/"+str(T_avg[0][0])+"--->"+str(T_avg[1][2])+"/"+str(T_avg[1][1])+"/"+str(T_avg[1][0]))

# Input name
filename=['simu_ctrl0_1d_20160101_grid_T']
print("Infile preamble: "+filename[0])

# Var info
varname=['sossheig']
varlongname=['sea_surface_height_above_geoid']
units=['[m]']
print("Var(s): "+varname[0]+" "+varlongname[0]+" "+units[0])

# Loop on input files
for n in range(0, len(filename)):
# Experiment N. 1
      # Path/name.extension 
      file1=pathdir+filename[n]+'.nc'
      # echo input file
      print('File.nc:'+file1)
      # open input file as a netcdf in order to read time range and var dims
      if os.path.isfile(file1)
        fh1=ncdf.Dataset(file1,'r')
        # select the var varname 
        var1 = fh1.variables[varname[n]][:]
        # read var dims
        aa=shape(var1)
        # select 3d var on sfc
        if len(aa)==3:
          VAR1=np.squeeze(var1[:,0,0])
        # select 4d var on sfc 
        if len(aa)==4:
          VAR1=np.squeeze(var1[:,0,0,0])
        # read time values
        time1 = fh1.variables['time'][:]
        # close the nc file
        fh1.close()
      else
        print ('File NOT found...Why?')
#Evaluates the time period
        # Compute days num = num of elements in time1 
        ndays=len(time1)
        print(time1[0])
        print(time1[-1])
        #
        # Defn of ini and end date from time var..
        data_ini=datetime.date(1950, 1, 1) + datetime.timedelta(days=time1[0]/24)
        data_fin=datetime.date(1950, 1, 1) + datetime.timedelta(days=time1[-1]/24+1)
#        #
#        times1 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")
## Evaluates the average value inthe specified period (T_avg)        
#        idx_avg_ini=(datetime.date(T_avg[0,0],T_avg[0,1],T_avg[0,2]) - data_ini).days
#        idx_avg_fin=(datetime.date(T_avg[1,0],T_avg[1,1],T_avg[1,2]) - data_ini).days
#        MEAN1=np.mean(VAR1[idx_avg_ini:idx_avg_fin])
#
## Experiment N. 2
#        file2=pathdir2+filename[n]+'.nc'
#        print(file2)
#        fh2=ncdf.Dataset(file2,'r')
#        var2 = fh2.variables[varname[n]][:]
#        aa=shape(var2)
#        if len(aa)==3:
#          VAR2=np.squeeze(var2[:,0,0])
#        if len(aa)==4:
#          VAR2=np.squeeze(var2[:,0,0,0])
#        time2 = fh2.variables['time'][:]
#        fh2.close()
##Evaluates the time period
#        ndays=len(time2)
#        data_ini=datetime.date(1900, 1, 1) + datetime.timedelta(days=time2[0]/24)
#        data_fin=datetime.date(1900, 1, 1) + datetime.timedelta(days=time2[-1]/24+1)
#        times2 = pd.date_range(str(data_ini),periods=ndays, freq = "1d")
## Evaluates the average value inthe specified period (T_avg)
#        idx_avg_ini=(datetime.date(T_avg[0,0],T_avg[0,1],T_avg[0,2]) - data_ini).days
#        idx_avg_fin=(datetime.date(T_avg[1,0],T_avg[1,1],T_avg[1,2]) - data_ini).days
#        MEAN2=np.mean(VAR2[idx_avg_ini:idx_avg_fin])
#
##PLOTS
#        fig = plt.figure(figsize=(9,5.5))
#        ax = fig.add_subplot(111)
#        ax.plot(times1,VAR1,'-r',label=expname[0]+': AVG= '+str(MEAN1)+units[n],linewidth=1.5)
#        ax.plot(times2,VAR2,'-b',label=expname[1]+': AVG= '+str(MEAN2)+units[n],linewidth=1.5)
#        leg = plt.legend(loc='upper left', ncol=2,  shadow=True, fancybox=True, fontsize=12)
#        leg.get_frame().set_alpha(0.3)
#        ax.grid('on')
#        plt.title( varlongname[n]+' '+units[n] ,fontsize=18)
#        ax.xaxis.set_major_locator(mdates.YearLocator())
#        ax.xaxis.set_minor_locator(mdates.MonthLocator((1,4,7,10)))
#        ax.xaxis.set_major_formatter(mdates.DateFormatter("\n%Y"))
#        ax.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
#        plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
#        plt.savefig(figdir+filename[n]+expname[0]+'_'+expname[1]+'.jpg')
#

