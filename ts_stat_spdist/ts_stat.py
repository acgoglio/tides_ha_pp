#
# Script for STATISTICS/SPECTRA/DISTRIB analisys on pnt time series
# You are supposed to provide a single file with the time series for each location. And the time-series is supposed to be obtained subtracting the mean value on the choosen period. This can be done by using the pextr.sh shell script.  
#
# by AC Goglio October 2019
#
# imports
import os # File and directory handling
import time # for waiting n seconds
#
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter 
#
from numpy import ma
from scipy.stats import ks_2samp # Kolmogorov-smirnof test
from statsmodels.distributions.empirical_distribution import ECDF # empirical distribution functions
#import rolling
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run paraeters:
# work dir path (WARNING: file in the directory will be removed if removing code is activated..)
workdir_path = '/work/ag15419/tmp/filters_Monaco/'
#
# dates
inidate = '01/07/2017'
enddate = '31/12/2017'
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Time interval to be used: ',dates_label)
# field(s) to be analized (to be improved with an array of fields..)
field_name='sossheig'
#
#
# OUTPUTS
# Flags for analisys type and output names
#------------------
# for ts statistics
stat_flag=0
stat_outfile_pre='stat_'
#---------------------------------
# for spectra computation and plot
spect_flag=1
spect_outfile_pre='spt_OBS_T8_f'
#spect_outfile_pre='spt_ctrl_'
spect_log=1 # Log y scale + top scale with period in hours (=1) or not (=0)
spt_color='magenta' # MOD color 
spt_obscolor='navy' # OBS color or first model color (for 2 model spetctra)
#
rate_s = 3600 # rate of data collection in seconds (3600 for hourly datasets; 86400 for daily datasets )
# Reduced range flag (venezia acqua alta event)
redrange_flag=0
# Smoothing spectra flag
smooth_flag=1
#-----------------------------
# For error (MOD-OBS) distribution comparison 
#er_distrib_flag=0
#er_distrib_outfile_pre='er_dst_'
#
# For comparison between time series values
val_distrib_flag=0
val_distrib_outfile_pre='val_t8ctrl_'

#---------------------------
print('You choose the following RUN FLAGS:','stat_flag (for statistics)=',stat_flag,' spect_flag (for spectra)=',spect_flag,' val_distrib_flag (for bias and rmseof (MOD-OBS) distributions)=',val_distrib_flag)
#
#
# INPUTS
# Path and name of inputs datasets
# Currently the extraction of ts in pnts is done externally by another script
#
# Flag for available dataset selection:
#
# OBSERVATION DATASET
# (csv and nc should contain the same ts in different locations)
csv_flag=0
nc_flag=1
# CSV obs file template: %csvobs_path%/%csvobs_prename%_stn_%csvobs_postname%.csv
# stn stands for the name of the location (can be a number, a name, etc.) the lists are below
if csv_flag == 1:
   print ('CSV obs dataset..')
   csvobs_path=workdir_path
   csvobs_fileprename='obs_'
   csvobs_postname=''
   csvobs_label='OBS ISPRA TG'
   #csv_stations_obs = ['ancona','carloforte','catania','crotone','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
   #csv_stations_lab = ['Ancona','Carloforte','Catania','Crotone','Imperia','Lampedusa','Livorno','Messina','Ortona','Palinuro','Trieste','Venezia','Vieste']
   #csv_stations_obs = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
   #csv_stations_lab = ['Ancona','Carloforte','Catania','Imperia','Lampedusa','Livorno','Messina','Ortona','Palinuro','Trieste','Venezia','Vieste']
   csv_stations_obs = ['napoli','civitavecchia','venezia','trieste'] 
   csv_stations_lab = ['napoli','civitavecchia','venezia','trieste']
   #
   print ('CSV obs file template: ',csvobs_path,csvobs_fileprename,'%stn_name%',csvobs_postname,'.csv')
   print('CSV stns: ',csv_stations_lab)
else:
   print ('No obs.csv dataset..')
   csv_stations_obs=[]
   csv_stations_lab=[]

# NC obs file template: %ncobs_path%/%ncobs_prename%_stn_%ncobs_postname%.nc
if nc_flag == 1:
   print ('NC obs dataset..')
   ncobs_path=workdir_path
   ncobs_fileprename=''
   ncobs_postname='TG_obs_HDF' #'TG_obs_okok_i' , TG_obs, TG_obsDF_i
   ncobs_label='Detided OBS' #EMODnet TG
   #nc_stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
   #nc_stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
   nc_stations_obs = ['Monaco','napoli','civitavecchia','Toulon','Tarifa']
   nc_stations_lab = ['Monaco','napoli','civitavecchia','Toulon','Tarifa']
   #
   print ('NC obs file template: ',ncobs_path,ncobs_fileprename,'%stn_name%',ncobs_postname,'.nc')
   print('NC stns: ',nc_stations_lab)
   
   if stat_flag == 1:
      diff_modnc_path = workdir_path
      diff_modnc_fileprename = 'diff_modnc_'
      diff_modnc_postnames = [ '_Tides8', '_Control_run' ]
      # Define arrays name dinamically
      for model_run in range(0,len(diff_modnc_postnames)):
          bias_array_name='bias_modnc_'+str(model_run)
          globals()[bias_array_name]=[]
          rmse_array_name='rmse_modnc_'+str(model_run)
          globals()[rmse_array_name]=[]
          # Num of obs per stn
          br_array_obsnum_name='br_obsnum_'+str(model_run)
          globals()[br_array_obsnum_name]=[]
else:
   print ('No obs.nc dataset..')
   nc_stations_obs=['trieste']
   nc_stations_lab=['trieste']
   diff_modnc_postnames=['trieste']
#
# MODEL DATASET
# Currently this scripts compute stats etc for obs vs model (just 1 dataset)
# Model 1st db file template (at the moment just one model db is possible): %model1obs_path%/%model1obs_prename%_stn_%model1obs_postname%.nc
num_of_models = 1
if num_of_models == 1:
   print ('Model dataset..')
   model_path=workdir_path
   model_fileprename=''
   model_postname='_mod_HDF_run' #_mod_Tides8 _mod_HDetided_run _mod_Tides8 _mod_Control_run
   model_label='Detided mod' #Tidal run
   #model_postname='_mod_Control_run'
   #model_label='Control run'
   #
   print ('Model file template: ',model_path,model_fileprename,'%stn_name%',model_postname,'.nc')

elif num_of_models == 2:
   # WARNING: in order to produce spectra of 2 model datasets ncobs_flag must be activated
   print ('Model datasets..')
   model_path=workdir_path
   model_fileprename=['','']
   model_postname=['_mod_Tides8','_mod_Tides8_original'] # '_mod_Tides8_run2yr','_mod_Detided_run2yr'
   model_label=['MOD (nearest interp)','MOD (4 points interp)'] # 'Tidal run','Detided run'
   model_labels=model_label
   rate_s_2 = rate_s # By default rate_s_2 = rate_s because data rates is the same in both datasets; for hourly rate_s_2 = 3600

elif num_of_models == 3:
   print ('Model datasets..')
   model_path=workdir_path
   model_fileprename=['','','']
   model_postname=['_mod_Tides8h','_mod_Detided_run','_mod_Control_run']
   model_label=['Tides_8 run','Detided run','Control run']

elif num_of_models == 4:
   print ('Model datasets..')
   model_path=workdir_path
   model_fileprename=['','','']
   model_postname=['_mod_Tides8','_mod_Control_run','_mod_HDetided_run']
   model_label=['Tides_8 run','Detided run','Control run']



# Required tidal constants
# Period defn
tidal_const_labels=['M2','K1','O1','S2','P1','N2','Q1','K2']
tidal_const_freqs=[]
# M2
M2freq = 1.0 / ( 12.421 * rate_s )
tidal_const_freqs.append(M2freq)
#K1
K1freq = 1.0 / (23.934 * rate_s )
tidal_const_freqs.append(K1freq)
#O1
O1freq = 1.0 / (25.819 * rate_s )
tidal_const_freqs.append(O1freq)
#S2
S2freq = 1.0 / (12.0 * rate_s )
tidal_const_freqs.append(S2freq)
#P1
P1freq = 1.0 / (24.066 * rate_s )
tidal_const_freqs.append(P1freq)
#N2
N2freq = 1.0 / (12.658 * rate_s )
tidal_const_freqs.append(N2freq)
#Q1
Q1freq = 1.0 / (26.868 * rate_s )
tidal_const_freqs.append(Q1freq)
#K2
K2freq = 1.0 / (11.967 * rate_s )
tidal_const_freqs.append(K2freq)

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################

# Move to the work dir and clean it if it exists otherwise create it
if os.path.exists(workdir_path) : 
   os.chdir(workdir_path) # mv to work dir
   print('I am moving in the directory: ',os.getcwd()) # pwd
   os.listdir() # ls
   # We do not want to remove input files if externally produced!!!!
   #print('WARNING: I am going to clean this directory in a while..')
   #time.sleep(10) # sleep for 10 seconds before removig everything in the work dir!
   #file2berm=os.path.join(workdir_path,'*.*') 
   #os.rm(file2berm) # clean dir
   #print ('Cleaning:',file2berm)
   #print('I am in the clean directory: ',os.getcwd()) # pwd
else:
  os.mkdir(workdir_path)
  os.chdir(workdir_path)
  print('I am in a new work directory: ',os.getcwd()) 
  

# Merge stations dbs and define models stz
print ('Merging the obs locations for model evaluation..')
# Merge all the obs stations
stations_obs = np.append(nc_stations_obs,csv_stations_obs)
# Define model stations (taking into account the "TG" string in EMODnet file names!)
nc_stations_obs_TG = [ str(nc_stations_obs_el)+'TG' for nc_stations_obs_el in nc_stations_obs ] 
#print ('nc_stations_obs_TG',nc_stations_obs_TG)
stations_mod = np.append(nc_stations_obs_TG,csv_stations_obs) 
obs_stations_lab = np.append(nc_stations_lab,csv_stations_lab)
mod_stations_lab = obs_stations_lab
print ('Model stns: ',mod_stations_lab)

# stn num
nc_numsta_obs=len(nc_stations_obs)
csv_numsta_obs=len(csv_stations_obs)
numsta_obs=len(stations_obs)
numsta_mod=len(stations_mod)
print ('obs stns # (# nc + # csv) = ',nc_numsta_obs,'+',csv_numsta_obs)
print ('mod stns # = ',numsta_mod)

# Check on stn numbers
if numsta_mod != numsta_obs:
   print ('WARNING: different number of stations in model/obs datasets! Why?!')

##################################
# Reading loop on stns of OBS dbs 
#################################

# ===============================
# CSV OBSERVATIONS DATASETS READING :
# ===============================
if csv_flag == 1:
    print ('Working on obs.csv dataset...')
    last_stn_idx=len(csv_stations_obs)
    print ('# of obs.csv stns = ',last_stn_idx)
    for stn in range(0,last_stn_idx):
        print('Reading values of CSV STN:',csv_stations_lab[stn])
        # Read values from CSV csv
        csv_pathname = csvobs_path+csvobs_fileprename+csv_stations_obs[stn]+'.csv'
        print ('csv obs file: ',csv_pathname )
        csv_obs = pd.read_csv(csv_pathname,sep=';',usecols=['idNum', 'station', 'year', 'month', 'day', 'hour', 'minu', 'sec', 'value'])
        print('...Done!')
        # Store csv obs values  
        csv_obs_values = csv_obs.value[:]
        # csv obs time in sec
        csv_te_obs= csv_obs_values.shape[0]
        csv_ts_obs= 0
        csv_time_obs=[]
        obscsv_val=[]
        mod4csv_val=[] # Array extracted from the model with the same elements of not-nan obs array 
        num_of_nan=0
        for csv_otime in range (csv_ts_obs,csv_te_obs):
            #csv_time_obs_el = (datetime.datetime(csv_obs.year[csv_otime],csv_obs.month[csv_otime],csv_obs.day[csv_otime],csv_obs.hour[csv_otime],csv_obs.minu[csv_otime])-datetime.datetime(inidate[6:10], inidate[3:5], inidate[0:2], 0, 0, 0)).total_seconds()/3600
            csv_time_obs_el = (datetime.datetime(int(csv_obs.year[csv_otime]),int(csv_obs.month[csv_otime]),int(csv_obs.day[csv_otime]),int(csv_obs.hour[csv_otime]),int(csv_obs.minu[csv_otime]))-datetime.datetime(2016, 10, 15, 0, 0, 0)).total_seconds()/rate_s
            # 
            # Handling of nan values in csv dataset
            if str(csv_obs_values[csv_otime]) != 'nan' : # rm nan values from arrays
               csv_time_obs.append(csv_time_obs_el)
               obscsv_val.append(float(csv_obs_values[csv_otime]))
            else:
               num_of_nan=num_of_nan+1
               print ('# of nan in csv dbs = ',num_of_nan)

        # rm mean value from csv dataset (for model and nc obs dataset this operation is done by the extractor script)
        mean_csv=np.mean(obscsv_val)
        print ('mean val to be subtracted from ts:', mean_csv)
        obscsv_val=obscsv_val[:]-mean_csv
        print ('Tini-Tend csv obs db =  ',csv_ts_obs,csv_te_obs)
        # ===============================
        # MODEL DATASET READING :
        # ===============================
        if num_of_models == 1:
           print ('Working on MODEL dataset...')
           # Read values from the model files
           stnmod=nc_numsta_obs+stn
           mod_pathname = model_path+model_fileprename+stations_mod[stnmod]+model_postname+'.nc'
           model = NC.Dataset(mod_pathname,'r')
           print ('model file: ',mod_pathname )
           # Store field values and do specific conversions if needed 
           if field_name == 'sossheig':
              print ('Working on field(s): ', field_name )
              print ('Conversion from meters to cm..')
              mod = model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
           # mod time in sec
           # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
           mod_time_offset = 852912.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
           time_mod = (model.variables["time_counter"][:]/rate_s)-mod_time_offset
           te_mod = mod.shape[0]
           ts_mod = 0
           mod_mean=np.mean(mod)
           print('Meants value to be subtracted: ',mod_mean)
           mod=mod[:]-mod_mean
           print ('Tini-Tend model dataset =  ',ts_mod,te_mod)

        elif num_of_models == 2:
           print ('Working on MODEL datasets...')
           mod_0=[]
           mod_1=[]
           for mox in 0,1:
              # Read values from the model files
              stnmod=nc_numsta_obs+stn
              mod_pathname=model_path+model_fileprename[mox]+stations_mod[stnmod]+model_postname[mox]+'.nc'
              model=NC.Dataset(mod_pathname,'r')
              print ('model file: ',mod_pathname )
              # Store field values and do specific conversions if needed 
              if field_name == 'sossheig':
                 print ('Working on field(s): ', field_name )
                 print ('Conversion from meters to cm..')
                 mod=model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
              # mod time in sec
              # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
              mod_time_offset = 591672.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0; for ??(acqua alta event in Venezia) = 852912.0 
              time_mod=(model.variables["time_counter"][:]/rate_s)-mod_time_offset
              te_mod=mod.shape[0]
              ts_mod=0
              mod_mean=np.mean(mod)
              print('Meants value to be subtracted: ',mod_mean)
              mod=mod[:]-mod_mean
              print ('Tini-Tend model dataset =  ',ts_mod,te_mod)
              str_mod='mod_'+str(mox)
              globals()[str_mod]=mod 

        elif num_of_models == 3:
           print ('Working on MODEL datasets...')
           mod_0=[]
           mod_1=[]
           mod_2=[]
           for mox in 0,1,2:
              # Read values from the model files
              stnmod=nc_numsta_obs+stn
              mod_pathname=model_path+model_fileprename[mox]+stations_mod[stnmod]+model_postname[mox]+'.nc'
              model=NC.Dataset(mod_pathname,'r')
              print ('model file: ',mod_pathname )
              # Store field values and do specific conversions if needed 
              if field_name == 'sossheig':
                 print ('Working on field(s): ', field_name )
                 print ('Conversion from meters to cm..')
                 mod=model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
              # mod time in sec
              # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
              mod_time_offset = 852912.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
              time_mod=(model.variables["time_counter"][:]/rate_s)-mod_time_offset
              te_mod=mod.shape[0]
              ts_mod=0
              mod_mean=np.mean(mod)
              print('Meants value to be subtracted: ',mod_mean)
              mod=mod[:]-mod_mean
              print ('Tini-Tend model dataset =  ',ts_mod,te_mod)
              str_mod='mod_'+str(mox)
              globals()[str_mod]=mod

        #######################################
        # TIME SERIES STATISTICS
        #######################################
        #if stat_flag == 1:

        #######################################
        # TIME SERIES SPECTRA
        #######################################
        if spect_flag == 1:

           spt_csv_obs=[]
           freq_csv_obs=[]
           freq_mod=[]
           spt_mod=[]
           #
           if num_of_models == 2:
              spt_csv_obs=abs(np.fft.fft(mod_0))
              freq_csv_obs = abs(np.fft.fftfreq(mod_0.shape[-1],rate_s))
              spt_mod=abs(np.fft.fft(mod_1))
              freq_mod = abs(np.fft.fftfreq(mod_1.shape[-1],rate_s))       
              csvobs_label=model_label[0]
              model_label=model_label[1]

           else:
              spt_csv_obs=abs(np.fft.fft(obscsv_val))
              freq_csv_obs = abs(np.fft.fftfreq(len(obscsv_val),rate_s))
              spt_mod=abs(np.fft.fft(mod))
              freq_mod = abs(np.fft.fftfreq(mod.shape[-1],rate_s))

           # Spectrum Plots for ISPRA dataset
           plt.figure(figsize=(15,8))
           plt.rc('font', size=14)
           plt.grid ()
           if redrange_flag != 1:
              plt.xlim(0.000005,0.00004)
           else:
              plt.xlim(0.000005,0.000025)
           plt.xlabel ('Frequency [Hz]')
           if spect_log == 1 :
              plt.yscale('log')
              plt.ylabel ('Log (Spectrum Amplitude)')
              plt.ylim(1,100000)
           elif spect_log == 0 :
              plt.ylabel ('Spectrum Amplitude')
              plt.ylim(1,60000)
           plt.title ('SSH'+' Spectrum - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label,pad=20)
           if redrange_flag != 1:
              # Add tidal constituents freqs
              if spect_log == 0 :
                 text_vertical_position=[50000,48000,50000,48000,44000,42000,42000,44000]
              elif spect_log == 1 :
                   text_vertical_position=[10,8,10,8,4,2,2,4]
              vlines_colors=['black','black','green','green','green','green','black','black']
              for idxx in range(0,8):
               plt.axvline(x=tidal_const_freqs[idxx], color=vlines_colors[idxx], linestyle="dashed")
               plt.text(tidal_const_freqs[idxx],text_vertical_position[idxx],tidal_const_labels[idxx] , color=vlines_colors[idxx] ,rotation=0,size=14)
               #
           else:
              if spect_log == 0 :
                 text_vertical_position=[50000,48000,50000,48000,44000,42000,42000,44000]
              elif spect_log == 1 :
                   text_vertical_position=[10,8,10,8,4,2,2,4]
              vlines_colors=['black','black','green','green','green','green','black','black']
              for idxx in range(0,8):
               plt.axvline(x=tidal_const_freqs[idxx], color=vlines_colors[idxx], linestyle="dashed")
               plt.text(tidal_const_freqs[idxx],text_vertical_position[idxx],tidal_const_labels[idxx] , color=vlines_colors[idxx] ,rotation=0,size=14)
           #
           plt.plot(freq_csv_obs,spt_csv_obs,color=spt_obscolor,label=csvobs_label)
           plt.plot(freq_mod,spt_mod,color=spt_color,label=model_label)
           # Add hour axes 
           if redrange_flag != 1:
              freqs2hours_labels=['55.5 h','27.8 h','18.5 h','13.9 h','11.1 h','9.3 h','7.9 h','6.9 h']
              freqs2hours_xposition=[0.000005,0.000010,0.000015,0.000020,0.000025,0.000030,0.000035,0.000040]
              freqs2hours_yposition=[100000,100000,100000,100000,100000,100000,100000,100000]
              for hidx in range(0,8):
                  plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=14)
           else:
              #freqs2hours_labels=['55.5 h','27.8 h','24 h','18.5 h','13.9 h','12 h','11.1 h']
              freqs2hours_labels=['55 h','37 h','28 h','24 h','19 h','16 h','14 h','12 h','11 h']
              freqs2hours_xposition=[0.000005,0.0000075,0.000010,0.000011574,0.000015,0.0000175,0.000020,0.000023148,0.000025]
              freqs2hours_xposition=np.array(freqs2hours_xposition)-0.0000002
              freqs2hours_yposition=105500 #[101500,101500,101500,101500,101500,101500,101500]
              for hidx in range(0,9):
                  plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition,freqs2hours_labels[hidx],size=14)
           #
           plt.legend( loc='upper left' )
           if spect_log == 1 :
              plt.savefig(workdir_path+spect_outfile_pre+'log_'+stations_mod[stnmod]+'_'+dates_label+'.jpg')
           elif spect_log == 0 :
              plt.savefig(workdir_path+spect_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
           plt.clf()


        ########################################
        # VALUE EMPIRICAL DISTRIBUTION FUNCTIONS
        ########################################
        if val_distrib_flag == 1:
           if num_of_models == 1:
              ecdf_csv_obs=[]
              ecdf_mod=[]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #plt.subplot(1,2,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              ecdf_csv_obs = ECDF(obscsv_val)
              ecdf_mod = ECDF(mod)
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              #plt.axvline(x=0.0, color='black')
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label)
              plt.plot(ecdf_csv_obs.x,ecdf_csv_obs.y,color='blue',label=csvobs_label) 
              plt.legend( loc='upper left' )
              # Kolmogorov-Smirnof test
              #plt.subplot(1,2,2)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()
              #scipy.stats.ks_2samp(mod,obscsv_val)

           elif num_of_models == 2:
              ecdf_csv_obs=[]
              ecdf_mod=[]
              ecdf_mod2=[]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #plt.subplot(1,2,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              ecdf_csv_obs = ECDF(obscsv_val)
              ecdf_mod = ECDF(mod_0)
              ecdf_mod2 = ECDF(mod_1)
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              #plt.axvline(x=0.0, color='black')
              max_dist2_value=round(np.max(ecdf_mod2.x),1)
              max_dist_value=round(np.max(ecdf_mod.x),1)
              max_dist_obs=round(np.max(ecdf_csv_obs.x),1)
              #min_dist2_value=round(np.min(ecdf_mod2.x),2)
              #min_dist_value=round(np.min(ecdf_mod.x),2)
              #min_dist_obs=round(np.min(ecdf_csv_obs.x),2)
              #
              plt.plot(ecdf_mod2.x,ecdf_mod2.y,color='black',label=model_label[1]+' (# = '+str(len(mod_1))+', max = '+str(max_dist2_value)+')')
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label[0]+' (# = '+str(len(mod_0))+', max = '+str(max_dist_value)+')')
              plt.plot(ecdf_csv_obs.x,ecdf_csv_obs.y,color='blue',label=csvobs_label+' (# = '+str(len(obscsv_val))+', max = '+str(max_dist_obs)+')')
              plt.legend( loc='upper left' )
              # Kolmogorov-Smirnof test
              #plt.subplot(1,2,2)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()
              #scipy.stats.ks_2samp(mod,obscsv_val)

           elif num_of_models == 3:
              ecdf_csv_obs=[]
              ecdf_mod=[]
              ecdf_mod2=[]
              ecdf_mod3=[]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #plt.subplot(1,2,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              ecdf_csv_obs = ECDF(obscsv_val)
              ecdf_mod = ECDF(mod_0)
              ecdf_mod2 = ECDF(mod_1)
              ecdf_mod3 = ECDF(mod_2)
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              #plt.axvline(x=0.0, color='black')
              plt.plot(ecdf_mod3.x,ecdf_mod3.y,color='black',label=model_label[2]+' (# = '+str(len(mod_2))+'; ssh range = ['+str(round(ecdf_mod3.x[1],2))+','+str(round(ecdf_mod3.x[-1],2))+']'+')')
              plt.plot(ecdf_mod2.x,ecdf_mod2.y,color='orange',label=model_label[1]+' (# = '+str(len(mod_1))+'; ssh range = ['+str(round(ecdf_mod2.x[1]))+','+str(round(ecdf_mod2.x[-1]))+']'+')')
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label[0]+' (# = '+str(len(mod_0))+'; ssh range = ['+str(round(ecdf_mod.x[1]))+','+str(round(ecdf_mod.x[-1]))+']'+')')
              plt.plot(ecdf_csv_obs.x,ecdf_csv_obs.y,color='blue',label=csvobs_label+' (# = '+str(len(obscsv_val))+'; ssh range = ['+str(round(ecdf_csv_obs.x[1]))+','+str(round(ecdf_csv_obs.x[-1]))+']'+')')
              plt.legend( loc='upper left' )
              # Kolmogorov-Smirnof test
              #plt.subplot(1,2,2)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()
              #scipy.stats.ks_2samp(mod,obscsv_val)
        


# ===============================
# NC OBSERVATIONS DATASETS READING :
# ===============================
if nc_flag == 1:
    print ('Working on obs.nc dataset...')
    last_stn_idx=len(nc_stations_obs)
    for stn in range(0,last_stn_idx):
        print('NC STN:',nc_stations_lab[stn])
        # Read values from NC csv
        ncobs_pathname = ncobs_path+ncobs_fileprename+nc_stations_obs[stn]+ncobs_postname+'.nc'
        print ('nc obs file: ',ncobs_pathname )
        nc_obs=[]
        nc_obs_values=[]
        nc_time_obs=[]
        nc_obs = NC.Dataset(ncobs_pathname,'r')
        # Store nc obs values and do specific conversions if needed 
        if field_name == 'sossheig':
           print ('Working on field(s): ', field_name )
           print ('Conversion from meters to cm..')
           nc_obs_values=nc_obs.variables[field_name][:,0]*100.0
        # nc obs time in sec
        #nc_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
        nc_time_offset = 591671.5 # for 1 Jul 2016 = 582911.5; for 1 Jul 2017 = 591671.5
        mod_time_offset = 591672.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
        print ( 'nc time offset (hours of the first day from 01/01/1950): ', nc_time_offset ,' Expected to be 586512.0' )
        if spect_flag != 1 and num_of_mod != 2:
         nc_time_obs = (nc_obs.variables["TIME"][:]*24)-nc_time_offset # want hours (not days) from the first day of analysis, not from 1/1/1950
        nc_te_obs = nc_obs_values.shape[0]
        nc_ts_obs = 0
        # subtract mean value from time series
        mean_nc_obs=np.mean(nc_obs_values)
        print('Meants value to be subtracted: ',mean_nc_obs)
        nc_obs_values=nc_obs_values[:]-mean_nc_obs
        print ('Tini-Tend nc obs db =  ',nc_ts_obs,nc_te_obs)

        # ===============================
        # MODEL DATASET READING :
        # ===============================
        if num_of_models == 1:
           print ('Working on MODEL dataset...')
           # Read values from the model files
           stnmod=stn
           mod_pathname = model_path+model_fileprename+stations_mod[stnmod]+model_postname+'.nc'
           model=[]
           mod=[]
           time_mod=[]
           model = NC.Dataset(mod_pathname,'r')
           print ('model file: ',mod_pathname )
           # Store field values and do specific conversions if needed 
           if field_name == 'sossheig':
              print ('Working on field(s): ', field_name )
              print ('Conversion from meters to cm..')
              mod = model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
              print ('Prova:',mod)
              print ('Prova:',model)
           # mod time in sec
           # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
           mod_time_offset = 591672.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
           time_mod = (model.variables["time_counter"][:]/rate_s)-mod_time_offset
           te_mod = mod.shape[0]
           ts_mod = 0
           #print ('Prova:',mod)
           mod_mean=np.mean(mod)
           print('Meants value to be subtracted: ',mod_mean)
           mod=mod[:]-mod_mean
           print ('Tini-Tend model dataset =  ',ts_mod,te_mod)

        elif num_of_models == 2:
           print ('Working on MODEL datasets...')
           mod_0=[]
           mod_1=[]
           for mox in 0,1:
               # Read values from the model files
               stnmod=stn
               mod_pathname = model_path+model_fileprename[mox]+stations_mod[stnmod]+model_postname[mox]+'.nc'
               model=[]
               mod=[]
               time_mod=[]
               model = NC.Dataset(mod_pathname,'r')
               print ('model file: ',mod_pathname )
               # Store field values and do specific conversions if needed 
               if field_name == 'sossheig':
                  print ('Working on field(s): ', field_name )
                  print ('Conversion from meters to cm..')
                  mod = model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
               # mod time in sec
               # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
               mod_time_offset = 591672.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
               time_mod = (model.variables["time_counter"][:]/rate_s)-mod_time_offset
               te_mod = mod.shape[0]
               ts_mod = 0
               mod_mean=np.mean(mod)
               print('Meants value to be subtracted: ',mod_mean)
               mod=mod[:]-mod_mean
               print ('Tini-Tend model dataset =  ',ts_mod,te_mod)
               str_mod='mod_'+str(mox)
               globals()[str_mod]=mod

        elif num_of_models == 3:
           print ('Working on MODEL datasets...')
           mod_0=[]
           mod_1=[]
           mod_2=[]
           for mox in 0,1,2:
               # Read values from the model files
               stnmod=stn
               mod_pathname = model_path+model_fileprename[mox]+stations_mod[stnmod]+model_postname[mox]+'.nc'
               model=[]
               mod=[]
               time_mod=[]
               model = NC.Dataset(mod_pathname,'r')
               print ('model file: ',mod_pathname )
               # Store field values and do specific conversions if needed 
               if field_name == 'sossheig':
                  print ('Working on field(s): ', field_name )
                  print ('Conversion from meters to cm..')
                  mod = model.variables["sossheig"][:,0,0]*100.0 # want cm not meters
               # mod time in sec
               # mod_time_offset = (datetime.datetime(datetime.datetime(inidate[6:10],inidate[3:5],inidate[0:2], 0, 0, 0)-datetime.datetime(datetime.datetime(1950,1,1, 0, 0, 0)).total_seconds()/3600
               mod_time_offset = 591672.0 # for 1 gen 2015 = 569784.0; for 1 Jul 2016 = 582912.0; for 1 Jul 2017 = 591672.0
               time_mod = (model.variables["time_counter"][:]/rate_s)-mod_time_offset
               te_mod = mod.shape[0]
               ts_mod = 0
               mod_mean=np.mean(mod)
               print('Meants value to be subtracted: ',mod_mean)
               mod=mod[:]-mod_mean
               print ('Tini-Tend model dataset =  ',ts_mod,te_mod)
               str_mod='mod_'+str(mox)
               globals()[str_mod]=mod


        #######################################
        # TIME SERIES STATISTICS
        #######################################
        if stat_flag == 1:
           print ('STAT...')
           # Compute the diffs mod - ons.nc if not done externally or read it from file
           diff_modnc_ready=[]
           diff_modnc=[]
           model_idx=0
           for model_run in diff_modnc_postnames:
               diff_modnc_pathname = diff_modnc_path+diff_modnc_fileprename+stations_mod[stnmod]+model_run+'.nc'
               diff_modnc_ready = NC.Dataset(diff_modnc_pathname,'r')
               # Store field values and do specific conversions if needed 
               if field_name == 'sossheig':
                  print ('Working on field(s): ', field_name )
                  print ('Conversion from meters to cm..')
                  diff_modnc = diff_modnc_ready.variables["sossheig"][:,0,0]*100.0 # want cm not meters
               # BIAS computation from diffs
               print ('Computation of bias..')
               bias_arr_string='bias_modnc_'+str(model_idx)
               globals()[bias_arr_string].append(np.mean(diff_modnc))
               # RMSE comutation from diffs
               print ('Computation of rmse..')
               rmse_arr_string='rmse_modnc_'+str(model_idx)
               globals()[rmse_arr_string].append(np.sqrt((diff_modnc ** 2).mean()))
               # Computation of the number of obs
               print ('Counting nuo of obs..')
               br_array_obsnum_string='br_obsnum_'+str(model_idx)
               globals()[br_array_obsnum_string].append(len(diff_modnc))
               #
               model_idx=model_idx+1

        #######################################
        # TIME SERIES SPECTRA
        #######################################
        if spect_flag == 1:

           spt_nc_obs=[]
           freq_nc_obs=[]
           freq_mod=[]
           spt_mod=[]
           # Spectra comparison between 2 model datasets:
           if num_of_models == 2:
              spt_nc_obs=abs(np.fft.fft(mod_0))
              freq_nc_obs = abs(np.fft.fftfreq(mod_0.shape[-1],rate_s))
              spt_mod=abs(np.fft.fft(mod_1))
              freq_mod = abs(np.fft.fftfreq(mod_1.shape[-1],rate_s_2))
              ncobs_label=model_labels[0]
              model_label=model_labels[1]
           else: # Comparison between model and nc obs spectra
              spt_nc_obs=abs(np.fft.fft(nc_obs_values))
              freq_nc_obs = abs(np.fft.fftfreq(nc_obs_values.shape[-1],rate_s))
              spt_mod=abs(np.fft.fft(mod))
              freq_mod = abs(np.fft.fftfreq(mod.shape[-1],rate_s))

           # Spectrum Plots for EMODnet stns dataset
           plt.figure(figsize=(15,8))
           plt.rc('font', size=12)
           plt.grid ()

           # Hourly data
           if rate_s == 3600:
              ##plt.xlim(0.000005,0.00004)
              ##plt.xlim(0.0000001,0.000001)
              ##plt.xlim(0.0000001,0.00015)
              plt.xlim(0.0000001,0.00015)
              #plt.xlim(0.00001,0.00005)
              plt.xscale('log')
              plt.xlabel ('Log(Frequency [Hz])')
              if spect_log == 1 :
                 plt.yscale('log')
                 plt.ylabel ('Log (Spectrum Amplitude)')
                 #plt.ylim(1,100000)
                 plt.ylim(1,1000000)
              elif spect_log == 0 :
                 plt.ylabel ('Spectrum Amplitude')
                 #plt.ylim(1,60000)
                 plt.ylim(1,20000)
              plt.title ('SSH'+' Spectrum - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label, pad=20)
              # Add tidal constituents freqs
              if spect_log == 0 :
                 text_vertical_position=[50000,48000,50000,48000,44000,42000,42000,44000]
              elif spect_log == 1 :
                 text_vertical_position=[10,8,10,8,4,2,2,4]
              vlines_colors=['black','black','green','green','green','green','black','black']
              ###for idxx in range(0,8):
              ###    plt.axvline(x=tidal_const_freqs[idxx], color=vlines_colors[idxx], linestyle="dashed")
              ###    plt.text(tidal_const_freqs[idxx],text_vertical_position[idxx],tidal_const_labels[idxx] , color=vlines_colors[idxx] ,rotation=0,size=14)
   
              # Add hour axes 
              #freqs2hours_labels=['55.5 h','27.8 h','18.5 h','13.9 h','11.1 h','9.3 h','7.9 h','6.9 h']
              freqs2hours_labels=['50','40','30','25','20','15','17','12','10','9','8','7','Period [Hours]']
              freqs2hours_val=[50.0,40.0,30.0,25.0,20.0,15.0,17.0,12.0,10.0,9.0,8.0,7.0,6.9]
              freqs2sec=np.asarray(freqs2hours_val)*3600
              freqs2hours_xposition=1.0/(freqs2sec)
              #freqs2hours_xposition=[0.000005,0.000010,0.000015,0.000020,0.000025,0.000030,0.000035,0.000040]
              ##freqs2hours_yposition=[100000,100000,100000,100000,100000,100000,100000,100000]
              freqs2hours_yposition=[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000]
              ###for hidx in range(0,len(freqs2hours_labels)):
              ###    plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=12)
              ###    plt.axvline(x=freqs2hours_xposition[hidx], color='black')
                # Add day axes 
              #freqs2hours_labels=['365d','100d','28d','15d','10d','2d','24h','12h','2h']
              #freqs2hours_xposition=[0.0000000317,0.0000001157,0.0000004134,0.0000007716,0.000001157,0.000005787,0.00001157,0.00002315,0.0001389]
              freqs2hours_labels=['365d','100d','28d','15d','10d','8d','6d','4d','2d','24h','12h','10h','8h','6h','4h','2h']
              freqs2hours_xposition=[0.0000000317,0.0000001157,0.0000004134,0.0000007716,0.000001157,0.000001447,0.000001929,0.000002894,0.000005787,0.00001157,0.00002315,0.00002778,0.00003472,0.00004630,0.00006944,0.0001389]
              freqs2hours_yposition=[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000]
              for hidx in range(0,len(freqs2hours_labels)):
                  plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=12)
                  plt.axvline(x=freqs2hours_xposition[hidx], color='black') #, linestyle="dashed"

           #
           # Daily data
           elif rate_s == 86400:
              if rate_s_2 == rate_s:
                 #plt.xlim(0.0,0.00003)
                 plt.xlim(-0.01,0.125)
                 #plt.xlim(-0.01,0.51)
                 freq_nc_obs=freq_nc_obs*86400
                 freq_mod=freq_mod*86400
                 plt.xlabel ('Frequency [1/Days]')
                 if spect_log == 1 :
                    plt.yscale('log')
                    plt.ylabel ('Log (Spectrum Amplitude)')
                    plt.ylim(0.5,10000)
                    # Add day axes 
                    freqs2hours_yposition=[12000,10000,10000,10000,10000,10000,10000,10000,10000,10000]
                 elif spect_log == 0 :
                    plt.ylabel ('Spectrum Amplitude')
                    plt.ylim(0.5,2000)
                    # Add day axes
                    freqs2hours_yposition=[2050,2000,2000,2000,2000,2000,2000,2000,2000,2000]
                 #
                 plt.title ('SSH'+' Spectrum - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label, pad=20)
                 # Add tidal constituents freqs
                 if spect_log == 0 :
                    text_vertical_position=[50000,48000,50000,48000,44000,42000,42000,44000]
                 elif spect_log == 1 :
                    text_vertical_position=[10,8,10,8,4,2,2,4]
                 vlines_colors=['black','black','green','green','green','green','black','black']
                 #for idxx in range(0,8):
                 #    plt.axvline(x=tidal_const_freqs[idxx], color=vlines_colors[idxx], linestyle="dashed")
                 #    plt.text(tidal_const_freqs[idxx],text_vertical_position[idxx],tidal_const_labels[idxx] , color=vlines_colors[idxx] ,rotation=0,size=14)
    
                 # Add day axes 
                 freqs2hours_labels=['365','100','28','18','15','10','5','3.3','2.5','2 [Days]']
                 freqs2hours_xposition=[0.00274,0.01,0.036,0.056,0.067,0.1,0.2,0.3,0.4,0.5]
                 #freqs2hours_yposition=[2000,2000,2000,2000,2000,2000,2000,2000,2000,2000]
                 #
                 for hidx in range(0,10):
                     plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=12)
                     plt.axvline(x=freqs2hours_xposition[hidx], color='black', linestyle="dashed")
              #
              # For comparison between hourly and daily spectra
              else: 
                 plt.xlabel ('Log (Frequency [1/s])')
                 plt.xscale('log')
                 plt.xlim(0.00000001,0.00015)
                 #plt.xlim(0.0000001,0.000001)
                 if spect_log == 1 :
                    plt.yscale('log')
                    plt.ylabel ('Log (Spectrum Amplitude)')
                    plt.ylim(0.5,1000000)
                 elif spect_log == 0 :
                    plt.ylabel ('Spectrum Amplitude')
                    plt.ylim(0,1600)
                 plt.title ('SSH'+' Spectrum - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label, pad=20)

                 # Add day axes 
                 freqs2hours_labels=['365d','100d','28d','15d','2d','24h','12h','2h']
                 freqs2hours_xposition=[0.0000000317,0.0000001157,0.0000004134,0.0000007716,0.000005787,0.00001157,0.00002315,0.0001389]
                 freqs2hours_yposition=[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000]
                 for hidx in range(0,8):
                     plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=12)
                     plt.axvline(x=freqs2hours_xposition[hidx], color='black', linestyle="dashed")
           #
           if smooth_flag == 0:
              # Plot mod and obs spectrum
              plt.plot(freq_mod,spt_mod,color=spt_color,label=model_label)
              plt.plot(freq_nc_obs,spt_nc_obs,color=spt_obscolor,label=ncobs_label)
              #plt.plot(freq_mod,spt_mod,color=spt_color,label=model_label)
              ###plt.plot(freq_mod,spt_mod,'ro',label=model_label)
           # Spectra Smoothing
           else: 
              #def fun_curves(t,A0,A1,A2): #,A2,A3,A4):
              #    return (A0 * np.exp(-A1 * t) + A2 ) #(A0+A1*t+A2*t**2+A3*t**3+A4*t**4)
              #fitted_mod, cov_mod = curve_fit(fun_curves,freq_mod,spt_mod)
              #plt.plot(freq_mod, fun_curves(freq_mod,*fitted_mod),color=spt_color,linestyle=":")           
              #fitted_obs, cov_obs = curve_fit(fun_curves,freq_nc_obs,spt_nc_obs)
              #plt.plot(freq_nc_obs, fun_curves(freq_nc_obs,*fitted_obs),color=spt_obscolor,linestyle=":")

              def moving_average(x, w):
                  return np.convolve(x, np.ones(w), 'valid') / w
              num_avgpoints=11
              num_avgint=(num_avgpoints-1)/2
              mod_avg=moving_average(spt_mod, num_avgpoints)
              nc_obs_avg=moving_average(spt_nc_obs, num_avgpoints)
              # For writing curves value
              #for iidx,jidx in zip(freq_mod,spt_mod):
              #    if jidx > 20.0:
              #       plt.annotate(str(jidx),xy=(iidx,jidx))
              # Curve integral:
              #print (np.sum(spt_mod))
              #
              plt.plot(freq_nc_obs[int(num_avgint):-int(num_avgint)],nc_obs_avg,color=spt_obscolor,label=ncobs_label)
              plt.plot(freq_mod[int(num_avgint):-int(num_avgint)],mod_avg,color=spt_color,label=model_label)

           plt.legend( loc='upper right' )
           if spect_log == 1 :
              plt.savefig(workdir_path+spect_outfile_pre+'log_'+stations_mod[stnmod]+'_'+dates_label+'.jpg')
           elif spect_log == 0 :
              plt.savefig(workdir_path+spect_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
           plt.clf()

	#########################################
        # VALUES EMPIRICAL DISTRIBUTION FUNCTIONS
        #########################################
        if val_distrib_flag == 1:

           if num_of_models == 1:
              ecdf_nc_obs=[]
              ecdf_mod=[]
              # EDFs computation
              ecdf_nc_obs = ECDF(nc_obs_values)
              ecdf_mod = ECDF(mod)
              # To be implemented a way of computing diffs between EDFs.. 
              #for xidx in ecdf_nc_obs.x:
              #    ecdf_diff=ecdf_mod[ecdf_nc_obs.x,:]-ecdf_nc_obs[ecdf_nc_obs.x,:]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #
              #plt.subplot(2,1,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              #plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label)
              plt.plot(ecdf_nc_obs.x,ecdf_nc_obs.y,color='blue',label=ncobs_label)
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label) 
              plt.legend( loc='upper left' )
              #
              ## Kolmogorov-Smirnof test
              #plt.subplot(2,1,2)
              #x_toplot=[]
              #y_toplot=[]
              ##idx_toplot=0
              #print (len(ecdf_nc_obs.x)-1)
              #print (len(ecdf_mod.x)-1)
              #for mod_ks_idx in range(0,len(ecdf_mod.x)):
              #    for nc_ks_idx in range(0,len(ecdf_nc_obs.x)):
              #        if ecdf_mod.x[mod_ks_idx] == ecdf_nc_obs.x[nc_ks_idx]:
              #           x_toplot.append(ecdf_nc_obs.x[nc_ks_idx])
              #           y_toplot.append(ecdf_mod.y[mod_ks_idx]-ecdf_nc_obs.y[nc_ks_idx])
              #       #idx_toplot=idx_toplot+1
              #plt.plot(x_toplot,y_toplot)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()

           elif num_of_models == 2:
              ecdf_nc_obs=[]
              ecdf_mod=[]
              ecdf_mod2=[]
              # EDFs computation
              ecdf_nc_obs = ECDF(nc_obs_values)
              ecdf_mod = ECDF(mod_0)
              ecdf_mod2 = ECDF(mod_1)
              # To be implemented a way of computing diffs between EDFs.. 
              #for xidx in ecdf_nc_obs.x:
              #    ecdf_diff=ecdf_mod[ecdf_nc_obs.x,:]-ecdf_nc_obs[ecdf_nc_obs.x,:]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #
              #plt.subplot(2,1,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              plt.plot(ecdf_mod2.x,ecdf_mod2.y,color='black',label=model_label[1]+' (# = '+str(len(mod_1))+')')
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label[0]+' (# = '+str(len(mod_0))+')')
              plt.plot(ecdf_nc_obs.x,ecdf_nc_obs.y,color='blue',label=ncobs_label+' (# = '+str(len(nc_obs_values))+')')
              plt.legend( loc='upper left' )
              #
              ## Kolmogorov-Smirnof test
              #plt.subplot(2,1,2)
              #x_toplot=[]
              #y_toplot=[]
              ##idx_toplot=0
              #print (len(ecdf_nc_obs.x)-1)
              #print (len(ecdf_mod.x)-1)
              #for mod_ks_idx in range(0,len(ecdf_mod.x)):
              #    for nc_ks_idx in range(0,len(ecdf_nc_obs.x)):
              #        if ecdf_mod.x[mod_ks_idx] == ecdf_nc_obs.x[nc_ks_idx]:
              #           x_toplot.append(ecdf_nc_obs.x[nc_ks_idx])
              #           y_toplot.append(ecdf_mod.y[mod_ks_idx]-ecdf_nc_obs.y[nc_ks_idx])
              #       #idx_toplot=idx_toplot+1
              #plt.plot(x_toplot,y_toplot)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()

           elif num_of_models == 3:
              ecdf_nc_obs=[]
              ecdf_mod=[]
              ecdf_mod2=[]
              # EDFs computation
              ecdf_nc_obs = ECDF(nc_obs_values)
              ecdf_mod = ECDF(mod_0)
              ecdf_mod2 = ECDF(mod_1)
              ecdf_mod3 = ECDF(mod_2)
              # To be implemented a way of computing diffs between EDFs.. 
              #for xidx in ecdf_nc_obs.x:
              #    ecdf_diff=ecdf_mod[ecdf_nc_obs.x,:]-ecdf_nc_obs[ecdf_nc_obs.x,:]
              #
              plt.figure(figsize=(14,8))
              plt.rc('font', size=12)
              plt.grid ()
              #
              #plt.subplot(2,1,1)
              plt.title ('Empirical Distribution Function - TG: '+mod_stations_lab[stnmod]+' - DT: '+dates_label)
              plt.xlabel ('sossheig [cm]')
              plt.ylabel ('EDF')
              plt.axhline(y=0.5, color='black', linestyle="dashed")
              plt.plot(ecdf_mod3.x,ecdf_mod3.y,color='black',label=model_label[2]+' (# = '+str(len(mod_2))+'; ssh range = ['+str(round(ecdf_mod3.x[1],2))+','+str(round(ecdf_mod3.x[-1],2))+']'+')')
              plt.plot(ecdf_mod2.x,ecdf_mod2.y,color='orange',label=model_label[1]+' (# = '+str(len(mod_1))+'; ssh range = ['+str(round(ecdf_mod2.x[1],2))+','+str(round(ecdf_mod2.x[-1],2))+']'+')')
              plt.plot(ecdf_mod.x,ecdf_mod.y,color='red',label=model_label[0]+' (# = '+str(len(mod_0))+'; ssh range = ['+str(round(ecdf_mod.x[1],2))+','+str(round(ecdf_mod.x[-1],2))+']'+')')
              plt.plot(ecdf_nc_obs.x,ecdf_nc_obs.y,color='blue',label=ncobs_label+' (# = '+str(len(nc_obs_values))+'; ssh range = ['+str(round(ecdf_nc_obs.x[1],2))+','+str(round(ecdf_nc_obs.x[-1],2))+']'+')')
              plt.legend( loc='upper left' )
              #
              ## Kolmogorov-Smirnof test
              #plt.subplot(2,1,2)
              #x_toplot=[]
              #y_toplot=[]
              ##idx_toplot=0
              #print (len(ecdf_nc_obs.x)-1)
              #print (len(ecdf_mod.x)-1)
              #for mod_ks_idx in range(0,len(ecdf_mod.x)):
              #    for nc_ks_idx in range(0,len(ecdf_nc_obs.x)):
              #        if ecdf_mod.x[mod_ks_idx] == ecdf_nc_obs.x[nc_ks_idx]:
              #           x_toplot.append(ecdf_nc_obs.x[nc_ks_idx])
              #           y_toplot.append(ecdf_mod.y[mod_ks_idx]-ecdf_nc_obs.y[nc_ks_idx])
              #       #idx_toplot=idx_toplot+1
              #plt.plot(x_toplot,y_toplot)
              #
              plt.savefig(workdir_path+val_distrib_outfile_pre+stations_mod[stnmod]+'_'+dates_label+'.jpg')
              plt.clf()


        #######################################
        # VALUES DISTRIBUTIONS
        #######################################
        #if val_distrib_flag == 1:

 
    # HERE the cycle on stns ends..(but we are still in emodnet obs flag!)

    if stat_flag == 1:
       ###############################
       # BIAS again.. (Sum and Plot )
       ###############################
   
       # Reduce the num of character in stns names for space reasons
       re_nc_stations_lab=nc_stations_lab
       lab=0
       for stz_nm in nc_stations_lab:
           re_nc_stations_lab[lab]=str(stz_nm[:3])
           lab=lab+1
       print ('re_nc_stations_lab',re_nc_stations_lab)
   
       # Plot for two model run dataset case:
       if len(diff_modnc_postnames) == 2 :
          print ('bias_modnc_0: ',bias_modnc_0)
          print ('bias_modnc_1: ',bias_modnc_1)
          print ('rmse_modnc_0: ',rmse_modnc_0)
          print ('rmse_modnc_1: ',rmse_modnc_1)
          print ('num of obs 0: ',br_obsnum_0)
          print ('num of obs 1: ',br_obsnum_1)
   
          # Diff between BIASes and RMSEs
          d_bias=np.asarray(bias_modnc_0)-np.asarray(bias_modnc_1)
          d_rmse=np.asarray(rmse_modnc_0)-np.asarray(rmse_modnc_1)
   
          # Bias and RMSE plots
          # BIAS
          plt.figure(figsize=(20,10))
          plt.rc('font', size=12)
          #
          plt.subplot(2,1,1)
          plt.title ('BIAS mod/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('BIAS SSH [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,bias_modnc_1,'g-o',label='run'+diff_modnc_postnames[1])
          plt.plot(re_nc_stations_lab,bias_modnc_0,'r-o',label='run'+diff_modnc_postnames[0])
          for idx_restn in range(0,len(re_nc_stations_lab)):
               val2write=str(round(bias_modnc_0[idx_restn],3))
               plt.text(re_nc_stations_lab[idx_restn],bias_modnc_0[idx_restn],val2write,color='black',size=14)
          plt.plot ([],[],' ',label='Num of obs = '+str(br_obsnum_0[0]))
          plt.legend( loc='upper left' )
          # bias diff
          plt.subplot(2,1,2)
          plt.title ('BIAS Tide_8/obs - BIAS Tide_4/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('BIAS DIFF [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,d_bias,color='blue',label='BIAS Tide8/obs - BIAS Tide4/obs ')
          plt.legend( loc='upper left' )
          plt.savefig(workdir_path+stat_outfile_pre+'bias_'+dates_label+'.jpg')
   
          # RMSE
          plt.figure(figsize=(20,10))
          plt.rc('font', size=12)
          #
          plt.subplot(2,1,1)
          plt.title ('RMSE mod/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('RMSE SSH [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,rmse_modnc_1,'g-o',label='run'+diff_modnc_postnames[1])
          plt.plot(re_nc_stations_lab,rmse_modnc_0,'r-o',label='run'+diff_modnc_postnames[0])
          plt.plot ([],[],' ',label='Num of obs = '+str(br_obsnum_0[0]))
          plt.legend( loc='upper left' )
          #
          plt.subplot(2,1,2)
          plt.title ('RMSE Tide_8/obs - RMSE Tide_4/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('RMSE DIFF [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,d_rmse,color='blue',label='RMSE Tide8/obs - RMSE Tide4/obs ')
          plt.legend( loc='upper left' )
          #
          plt.savefig(workdir_path+stat_outfile_pre+'rmse_'+dates_label+'.jpg')
          plt.clf()
   
       # Plot for 1 model run dataset case:
       if len(diff_modnc_postnames) == 1 :
          print ('bias_modnc_0: ',bias_modnc_0)
          print ('rmse_modnc_0: ',rmse_modnc_0)
          print ('num of obs 0: ',br_obsnum_0)   

          # Bias and RMSE plots
          # BIAS
          plt.figure(figsize=(20,10))
          plt.rc('font', size=12)
          #
          plt.subplot(2,1,1)
          plt.title ('BIAS mod/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('BIAS SSH [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,bias_modnc_0,'r-o',label='run'+diff_modnc_postnames[0])
          plt.plot ([],[],' ',label='Num of obs = '+str(br_obsnum_0[0]))
          plt.legend( loc='upper left' )
          # RMSE
          plt.subplot(2,1,2)
          plt.title ('RMSE mod/obs - VAR: '+field_name+' - DT: '+dates_label)
          plt.grid ()
          plt.ylabel ('RMSE SSH [cm]')
          plt.axhline(y=0, color='black')
          plt.plot(re_nc_stations_lab,rmse_modnc_0,'r-o',label='run'+diff_modnc_postnames[0])
          plt.plot ([],[],' ',label='Num of obs = '+str(br_obsnum_0[0]))
          plt.legend( loc='upper left' )
          #
          plt.savefig(workdir_path+stat_outfile_pre+dates_label+'.jpg')
          plt.clf()


