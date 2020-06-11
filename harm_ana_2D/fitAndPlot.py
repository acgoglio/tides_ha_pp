#
# Script for DETIDING
# by Salish Sea MEOPAR
# https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis/raw/tip/compare_tides/Analysis8Components.ipynb
#
# imports
#%matplotlib inline
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
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
grid = 'T' # Choose T, U, V or uv2t grid
num_of_models = 1 # 1 for dataset maps or 2 for diffs between 2 datasets
dataset_name = '8' # '8' for Tides_8 run, 'ctrl' for the Control run or 'diff' for computing the diffs
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
workdir_path = '/work/ag15419/tmp/HA_nn/ha_map/'+'/'
archive_path = '/work/ag15419/tmp/HA_nn/ha_map/'+'/'
path = '/work/ag15419/tmp/2d_harm_ana/'
# dates
#---------------------
# Choose start and end dates of the period 
inidate = '01/07/2017'
enddate = '31/12/2017'
#--------------------
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
#
# Fields to be analyzed
#--------------------
if grid == 'T':
   # 2D fields:
   var_2d='sossheig' # ['sossheig','sowaflup','soevapor','soprecip','sorunoff','soshfldo','sohefldo','solofldo','sosefldo','solafldo','somxl010']
   field_2d_units='m' #['m','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','W/m2','W/m2','W/m2','W/m2','m']
   vlev_val=0
   # Fix thereshold for max and min map values for 2d fields or not
   T_threshold_2d=0
   field_2d_inf=[10,36]
   field_2d_sup=[25,42]
#
# INPUTS
# Path and name of inputs datasets
# Currently the extraction of nc is done externally by another script
#
# MODEL DATASETS
# Currently this scripts plots maps for just 1 dataset (can be improoved..)
# Model 1st db file template (at the moment just one model db is possible): %model1obs_path%/%model1obs_prename%_stn_%model1obs_postname%.nc
if num_of_models == 1:
   if dataset_name == '8':
      model_path=archive_path
      model_fileprename='map'
      model_postname='mod_Tides8' #'mod_Tides8' or 'mod_Control_run'
      model_label='Tides_8 run' #'Tides_8 run' or 'Control run'
   elif dataset_name == 'ctrl':
      model_path=archive_path
      model_fileprename='map'
      model_postname='mod_Control_run' #'mod_Tides8' or 'mod_Control_run'
      model_label='Control run' #'Tides_8 run' or 'Control run'
   print ('Model file templates: ',model_path,model_fileprename,'?D_','%vlev%_','%field%',dates_label,'_',model_postname,'.nc')

# Outfile 
new_path=workdir_path
new_fileprename='amppha'


########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################
##
#Tidal outputs list (filds in the new nc file)
tides8_out=['M2_Amp','M2_Pha','K1_Amp','K1_Pha','O1_Amp','O1_Pha','S2_Amp','S2_Pha','P1_Amp','P1_Pha','N2_Amp','N2_Pha','Q1_Amp','Q1_Pha','K2_Amp','K2_Pha']
tides8_out_A=['M2_Amp','K1_Amp','O1_Amp','S2_Amp','P1_Amp','N2_Amp','Q1_Amp','K2_Amp']
tides8_out_P=['M2_Pha','K1_Pha','O1_Pha','S2_Pha','P1_Pha','N2_Pha','Q1_Pha','K2_Pha']
### Constants defn 

### M2
M2freq = 28.984106 # degrees per hour
M2freq = M2freq*np.pi/180. # radians per hour
###K1
K1freq = 15.041069*np.pi/180.
###O1
O1freq = 13.943036*np.pi/180.
###S2
S2freq = 30.000002*np.pi/180.
###P1
P1freq = 14.958932*np.pi/180.
###N2
N2freq = 28.439730*np.pi/180.
###Q1
Q1freq = 13.398661*np.pi/180.
###K2
K2freq = 30.082138*np.pi/180.
##

# Fitting
def octuple_mod(t, M2amp_mod, M2pha_mod, K1amp_mod, K1pha_mod, O1amp_mod, O1pha_mod, S2amp_mod, S2pha_mod, P1amp_mod, P1pha_mod, N2amp_mod, N2pha_mod, Q1amp_mod, Q1pha_mod, K2amp_mod, K2pha_mod):
    return (M2amp_mod*np.cos(M2freq*t-M2pha_mod*np.pi/180.)+
            K1amp_mod*np.cos(K1freq*t-K1pha_mod*np.pi/180.)+
            O1amp_mod*np.cos(O1freq*t-O1pha_mod*np.pi/180.)+
            S2amp_mod*np.cos(S2freq*t-S2pha_mod*np.pi/180.)+
            P1amp_mod*np.cos(P1freq*t-P1pha_mod*np.pi/180.)+
            N2amp_mod*np.cos(N2freq*t-N2pha_mod*np.pi/180.)+
            Q1amp_mod*np.cos(Q1freq*t-Q1pha_mod*np.pi/180.)+
            K2amp_mod*np.cos(K2freq*t-K2pha_mod*np.pi/180.))


# Open input files and read infos
nc2open=model_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'.nc'
# TMP:
#nc2open='/work/ag15419/tmp/2d_harm_ana_6m/amppha2D_0_sossheig_20170701_20171231_mod_Tides8.nc_ok'
print ('Input file = ',nc2open)
model = NC.Dataset(nc2open,'r')
nc2create=new_path+new_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'.nc'
amppha=NC.Dataset(nc2create,'a')
print ('Output file = ',nc2create)

# Read lat, lon and field values 
if grid == 'uv2t':
   lons = model.variables['lon'][:]
   lats = model.variables['lat'][:]
   time_mod=(model.variables['time_counter'][:]/3600.0)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
   vals = (model.variables[var_2d][:])
else:
   lons = model.variables['nav_lon'][:]
   lats = model.variables['nav_lat'][:]
   time_mod=(model.variables['time_counter'][:]/3600.0)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
   vals = (model.variables[var_2d][:])
   print (lats.shape,lons.shape)
# Sel vertical levs
print('vlev=',vlev_val)
# Initial and end time-steps
ts_mod=0
te_mod = vals.shape[0]

# Write lat lon and time infos in the new file
#lat_dim = nc2create.createDimension('lat', len(lats)) # latitude axis
#lon_dim = nc2create.createDimension('lon', len(lons)) # longitude axis
#time_dim = nc2create.createDimension('time', None) 
#nc2create.title='Amplitur=de and Phases'
#nc2create.subtitle="Obtained from a fit"
#lat = nc2create.createVariable('lat', np.float32, ('lat','lon'))
#lat.units = 'degrees_north'
#lat.long_name = 'latitude'
#lat=lats
#print(lats.shape, lats)
#lon = nc2create.createVariable('lon', np.float32, ('lat','lon'))
#lon.units = 'degrees_east'
#lon.long_name = 'longitude'
#lon=lons
#print(lons.shape, lons)
#time = nc2create.createVariable('time', np.float64, ('time'))
#time.units = 'hours since 1950-01-01'
#time.long_name = 'time'

# Create new fields:
###for idx_t8field in tides8_out_A:
###    globals()[idx_t8field] = nc2create.createVariable(globals()[idx_t8field],np.float64,('y','x')) 
###    idx_t8field.units = 'm' 
###    idx_t8field.standard_name = 'Amplitude' 
###for idx_t8field in tides8_out_P:
###    idx_t8field = nc2create.createVariable(idx_t8field,np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
###    idx_t8field.units = 'deg' # degrees 
###    idx_t8field.standard_name = 'Phase' # this is a CF standard name

M2_Amp=amppha.createVariable('M2_Amp',np.float64,('y','x'))
M2_Amp.units = 'm' 
M2_Amp.standard_name = 'Amplitude'
#
K1_Amp=amppha.createVariable('K1_Amp',np.float64,('y','x'))
K1_Amp.units = 'm' 
K1_Amp.standard_name = 'Amplitude'
#
O1_Amp=amppha.createVariable('O1_Amp',np.float64,('y','x'))
O1_Amp.units = 'm'
O1_Amp.standard_name = 'Amplitude'
#
S2_Amp=amppha.createVariable('S2_Amp',np.float64,('y','x'))
S2_Amp.units = 'm'
S2_Amp.standard_name = 'Amplitude'
#
P1_Amp=amppha.createVariable('P1_Amp',np.float64,('y','x'))
P1_Amp.units = 'm'
P1_Amp.standard_name = 'Amplitude'
#
N2_Amp=amppha.createVariable('N2_Amp',np.float64,('y','x'))
N2_Amp.units = 'm'
N2_Amp.standard_name = 'Amplitude'
#
Q1_Amp=amppha.createVariable('Q1_Amp',np.float64,('y','x'))
Q1_Amp.units = 'm'
Q1_Amp.standard_name = 'Amplitude'
#
K2_Amp=amppha.createVariable('K2_Amp',np.float64,('y','x'))
K2_Amp.units = 'm'
K2_Amp.standard_name = 'Amplitude'
#
#
#
M2_Pha=amppha.createVariable('M2_Pha',np.float64,('y','x'))
M2_Pha.units = 'deg'
M2_Pha.standard_name = 'Phase'
#
K1_Pha=amppha.createVariable('K1_Pha',np.float64,('y','x'))
K1_Pha.units = 'deg'
K1_Pha.standard_name = 'Phase'
##
O1_Pha=amppha.createVariable('O1_Pha',np.float64,('y','x'))
O1_Pha.units = 'deg'
O1_Pha.standard_name = 'Phase'
##
S2_Pha=amppha.createVariable('S2_Pha',np.float64,('y','x'))
S2_Pha.units = 'deg'
S2_Pha.standard_name = 'Phase'
##
P1_Pha=amppha.createVariable('P1_Pha',np.float64,('y','x'))
P1_Pha.units = 'deg'
P1_Pha.standard_name = 'Phase'
##
N2_Pha=amppha.createVariable('N2_Pha',np.float64,('y','x'))
N2_Pha.units = 'deg'
N2_Pha.standard_name = 'Phase'
##
Q1_Pha=amppha.createVariable('Q1_Pha',np.float64,('y','x'))
Q1_Pha.units = 'deg'
Q1_Pha.standard_name = 'Phase'
##
K2_Pha=amppha.createVariable('K2_Pha',np.float64,('y','x'))
K2_Pha.units = 'deg'
K2_Pha.standard_name = 'Phase'
#
print(vals.shape)
#Loop on grid points:
print ('len(lons): ',len(lons[1]),' len(lats): ',len(lats[0]))
for idx_lon in range (0,1307):  #len(lons[1])):
    for idx_lat in range (0,380):    #len(lats[0])):
 
        # Initialize the array storing Amp an Pha values
        fitted_mod=[]

        #longitude=lons[idx_lat,idx_lon]
        #latitude=lats[idx_lat,idx_lon]
        print ('Indexes and Coordinates (lon/lat):', idx_lon, idx_lat )
        #
        ssh_mod = vals[:,idx_lat,idx_lon]
        #print ('SSH time-series:',ssh_mod)
        #
        #print ("Fitting the model...")
        try:
           fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:te_mod],ssh_mod[ts_mod:te_mod])
           #print ('As and Ps:', fitted_mod)
           #print ('Amp M2: ', fitted_mod[0])
        except RuntimeError:
           print("Error - curve_fit failed")
           fitted_mod=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        # Write in the new nc file
        ###for idx_t8field in range(0,len(tides8_out)):
        ###    globals()[tides8_out[idx_t8field]][:,idx_lat,idx_lon]=fitted_mod[idx_t8field]
        M2_Amp[idx_lat,idx_lon]=fitted_mod[0]
        K1_Amp[idx_lat,idx_lon]=fitted_mod[2]
        O1_Amp[idx_lat,idx_lon]=fitted_mod[4]
        S2_Amp[idx_lat,idx_lon]=fitted_mod[6]
        P1_Amp[idx_lat,idx_lon]=fitted_mod[8]
        N2_Amp[idx_lat,idx_lon]=fitted_mod[10]
        Q1_Amp[idx_lat,idx_lon]=fitted_mod[12]
        K2_Amp[idx_lat,idx_lon]=fitted_mod[14]

        M2_Pha[idx_lat,idx_lon]=fitted_mod[1]
        K1_Pha[idx_lat,idx_lon]=fitted_mod[3]
        O1_Pha[idx_lat,idx_lon]=fitted_mod[5]
        S2_Pha[idx_lat,idx_lon]=fitted_mod[7]
        P1_Pha[idx_lat,idx_lon]=fitted_mod[9]
        N2_Pha[idx_lat,idx_lon]=fitted_mod[11]
        Q1_Pha[idx_lat,idx_lon]=fitted_mod[13]
        K2_Pha[idx_lat,idx_lon]=fitted_mod[15]

#print(amppha)
amppha.close()



