#
# Script for MAPS plot  
# You are supposed to provide a single nc file with the mean of the field on a single period and/or single file with field averaged on seasons. These can be built by map_extr.sh_nov shell script.  
#
# by AC Goglio November 2019
#
#source activate mappyenv # activate my virtual environment
#
# imports
import os # File and directory handling
import time # for waiting n seconds
#
import matplotlib as mpl # Palette
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
from mpl_toolkits.basemap import Basemap # For plotting maps
from matplotlib.colors import LogNorm # For log scale in pcolor 
import seaborn as sns # Darker Palette
from matplotlib.collections import LineCollection # Line between points
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
grid = 'uv2t' # Choose T, U, V or uv2t grid
num_of_models = 2 # 1 for dataset maps or 2 for diffs between 2 datasets
dataset_name = 'diff' # '8' for Tides_8 run, 'ctrl' for the Control run or 'diff' for computing the diffs, 'bathy'for bathymetry
sub_plot_flag = 0 # flag for sub area plots (if == 1 are produced the global and sub areas plots)
tidegauge_flag = 0 # To add locations and numbers of TG (in order to use this must be dataset_name = 'bathy')
transpsect_flag = 0 # To add sections where the transports are computed
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
#workdir_path = '/work/ag15419/tmp/map_ana_'+grid+dataset_name+'/'
workdir_path = '/work/ag15419/tmp/SSH_25_3/currents_extr/'
#workdir_path = '/work/ag15419/tmp/map_bathy/'
#
# dates
#---------------------
# Choose start and end dates of the period 
inidate = '01/01/2015'
enddate = '01/03/2015'
# Choose periods for the analysis (yr means on the whole period), single file with seasonal mean should exists!
period_of_analysis=['yr','DJF','MAM','JJA','SON'] # 'yr','DJF','MAM','JJA','SON','months'(to be implemented)
#--------------------
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
print ('Periods of analysis: ',period_of_analysis)
#
# field(s) to be plotted
#-----------------------
if grid == 'T':
   if dataset_name != 'bathy':
      # T grid
      mod_depth_var='deptht'
      # 3D fields:
      field_3d_name=['vosaline'] # #'votemper','vosaline'
      field_3d_units=['PSU'] # 'degC','PSU'
      # field_3d_lev_val is the index of the level , field_3d_lev is the depth in meters 
      field_3d_lev_val=[0,10,30] #[0,3,10,23,30,43,59,73,95] # Mod lev num
      field_3d_lev=[1,30,150] #[1,8,30,100,150,300,600,1000,2000] # Approx depth
      # Fix thereshold for max and min map values for 3d fields or not
      T_threshold_3d=0
      field_3d_inf=[10,36]
      field_3d_sup=[25,42]
      #
      Tdiff_threshold_3d=0
      tdiff_th=0.6
      # 2D fields:
      field_2d_name=[] #['sossheig','sowaflup','soevapor','soprecip','sorunoff','soshfldo','sohefldo','solofldo','sosefldo','solafldo','somxl010']
      field_2d_units=['m'] #['m','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','W/m2','W/m2','W/m2','W/m2','m']
      field_2d_lev=0
      # Fix thereshold for max and min map values for 2d fields or not
      T_threshold_2d=0
      field_2d_inf=[-10,10]
      field_2d_sup=[-10,10]
   else:
      # T grid
      mod_depth_var='deptht'
      # 3D fields:
      field_3d_name=[] # 
      field_3d_units=[] #
      # field_3d_lev_val is the index of the level , field_3d_lev is the depth in meters 
      field_3d_lev_val=[] #[0,3,10,23,30,43,59,73,95] # Mod lev num
      field_3d_lev=[] #[1,8,30,100,150,300,600,1000,2000] # Approx depth
      # Fix thereshold for max and min map values for 3d fields or not
      T_threshold_3d=0
      field_3d_inf=[]
      field_3d_sup=[]
      #
      Tdiff_threshold_3d=0
      tdiff_th=0.5
      # 2D fields:
      field_2d_name=['Bathymetry']
      field_2d_units=['m']
      field_2d_lev=0
      # Fix thereshold for max and min map values for 2d fields or not
      T_threshold_2d=1
      field_2d_inf=[-5000]
      field_2d_sup=[0]

elif grid == 'U':
     mod_depth_var='depthu'
     # U grid
     field_3d_name=['vozocrtx'] 
     field_3d_inf=[-0.4]
     field_3d_sup=[0.4]
     field_3d_units=['m/s']
     field_3d_lev_val=[0,3,10,23,30,43,59,73,95]
     field_3d_lev=[1,8,30,100,150,300,600,1000,2000]
     #
     Tdiff_threshold_3d=1
     tdiff_th=0.15
     #
     field_2d_name=['sozotaux']
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=['N/m2']
elif grid == 'V':
     mod_depth_var='depthv' 
     # V grid
     field_3d_name=[] # 'vomecrty'
     field_3d_inf=[-0.3]
     field_3d_sup=[0.3]
     field_3d_units=['m/s']
     field_3d_lev_val=[0,3,10,23,30,43,59,73,95]
     field_3d_lev=[1,8,30,100,150,300,600,1000,2000]
     #
     Tdiff_threshold_3d=1
     tdiff_th=0.2
     #
     field_2d_name=['sometauy']
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=['N/m2']
elif grid == 'uv2t':
     mod_depth_var='depth'
     # uv2t grid
     field_3d_name=['i']
     field_3d_inf=[0]
     field_3d_sup=[0.06]
     field_3d_units=['m/s']
     field_3d_lev_val=[4,35] #,10,23,30,43,59,73,95] #,59,73,95] #[0,3,10,23,30,43,59,73,95] or [4,35]
     field_3d_lev=[10,200] #30,100,150,300,600,1000,2000] #,600,1000,2000] #[1,8,30,100,150,300,600,1000,2000] or [10,200]
     #
     Tdiff_threshold_3d=1
     tdiff_th=0.01
     #
     field_2d_name=[]
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=[]
#
# SUB AREAS for currents analysis
# Gibraltar, Adriatic and Eastern Med Areas
Area_G_lon=[-7,-2.5] # GB LONs [-7,9]
Area_G_lat=[35,37] # GB LATs [34,45]
Area_A_lon=[11,23] # AD LONs 
Area_A_lat=[37,46] # AD LATs 
Area_EM_lon=[25.50,27.00] # EM LONs All:[9,37]  Dardanelles: 25.50,27.00
Area_EM_lat=[39.50,40.50] # EM LATs All: [9,37] Dardanelles: 39.50,40.50
Area_ME_lon=[14,17] # ME LONs
Area_ME_lat=[36.5,39] # ME LATs
#
Area_name_ar=['Gibraltar','Adriatic','Dardanelles','Messina'] #'Gibraltar','Adriatic','Eastern_Med','Messina'
Area_minlon_ar=[Area_G_lon[0],Area_A_lon[0],Area_EM_lon[0],Area_ME_lon[0]]
Area_maxlon_ar=[Area_G_lon[1],Area_A_lon[1],Area_EM_lon[1],Area_ME_lon[1]]
Area_minlat_ar=[Area_G_lat[0],Area_A_lat[0],Area_EM_lat[0],Area_ME_lat[0]]
Area_maxlat_ar=[Area_G_lat[1],Area_A_lat[1],Area_EM_lat[1],Area_ME_lat[1]]
#
# U and V ranges in each region
if grid == 'U':
   Area_vmin_ar=[-0.5,-0.2,-0.4,-0.5]
   Area_vmax_ar=[1.0,0.2,0.4,0.5]
elif grid == 'V':
   Area_vmin_ar=[-0.5,-0.3,-0.2,-0.5]
   Area_vmax_ar=[0.5,0.3,0.2,0.2]
elif grid == 'T' and dataset_name != 'bathy':
   Area_vmin_ar=[-0.5,-0.5,-0.8,2.0] # -1.0,-0.5,-0.8,2.0
   Area_vmax_ar=[0.5,0.5,0.8,-2.0] # 1.0,0.5,0.8,-2.0
elif grid == 'uv2t':
   Area_vmin_ar=[-0.05,-0.05,-0.05,-0.05]
   Area_vmax_ar=[0.05,0.05,0.05,0.05]
elif dataset_name == 'bathy':
   Area_vmin_ar=[-5000,-5000,-5000,-5000]
   Area_vmax_ar=[0.0,0.0,0.0,0.0]
# INPUTS
# Path and name of inputs datasets
# Currently the extraction of nc is done externally by another script
#
# MODEL DATASETS
# Currently this scripts plots maps for just 1 dataset (can be improoved..)
# Model 1st db file template (at the moment just one model db is possible): %model1obs_path%/%model1obs_prename%_stn_%model1obs_postname%.nc
if num_of_models == 1:
   if dataset_name == '8':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_Tides8' #'mod_Tides8' or 'mod_Control_run'
      model_label='Tides_8 run' #'Tides_8 run' or 'Control run'
   elif dataset_name == 'ctrl':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_Ctrl' #'mod_Tides8' or 'mod_Control_run'
      model_label='Control run' #'Tides_8 run' or 'Control run'
   elif dataset_name == 'bathy':
      model_path='/work/ag15419/PHYSW24_DATA/TIDES/DATA0/'
      model_fileprename='bathy'
      model_postname='meter' 
      model_label='Bathymetry'
   elif dataset_name == 'Tides8-25cm':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_Tides8-25cm'
      model_label='Tides_8-25cm' 
   elif dataset_name == 'Tides8-25VV':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_25VV'
      model_label='Simu_tides8_v6'
   elif dataset_name == 'Tides8-FES2014':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_FES2014'
      model_label='Simu_tides8_v7'
   elif dataset_name == 'Tides8_v8':
      model_path=workdir_path
      model_fileprename='map'
      model_postname='mod_Simu_tides8_v8'
      model_label='Simu_tides8_v8'
   # map3D_yr_%LEV%_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   if dataset_name != 'bathy':
      print ('Model file templates: ',model_path,model_fileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_',model_postname,'.nc')
   else:
      print ('Model file templates: ',model_path,model_fileprename,'_',model_postname,'.nc')
#
elif num_of_models == 2:
   print ('Model datasets DIFF MAPS..')
   model_path=workdir_path
   model_infileprename='map'
   model_fileprename='diff_map'
   model_postname=['mod_Simu_tides8_v8','mod_imu_tides8_v6'] #,'mod_Control_run'], 'mod_Tides8-25cm'] 
   model_label=['','( Simu_tides8_v8 - Simu_tides8_v6 )' ] #,'Control run'] 'Tides_8-25cm' ]
   #
   #model_postname=['mod_GBmod','mod_Control_run']
   #model_label=['GB mod','Control run']
   #
   #model_postname=['mod_GBmod','mod_Tides8']
   #model_label=['GB mod','Tides_8 run']
   #
   # map3D_yr_%LEV%_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   print ('Model file templates: ',model_path,model_infileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_','%model_postname%','.nc')

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################

# Move to the work dir and clean it if it exists otherwise create it
if os.path.exists(workdir_path) : 
   os.chdir(workdir_path) # mv to work dir
   print('I am moving in the directory: ',os.getcwd()) # pwd
   #os.listdir() # ls
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
  

# List input/output cases  
print ('Model dataset: ',model_label)
print ('Fields (3D+2D): ',field_3d_name,field_2d_name)
print ('# of 3d fields vertical levels: ',field_3d_lev)
print ('Periods (global+season): ',dates_label,period_of_analysis)

if num_of_models == 1:

    print ('I am working on the following model dataset: ',model_label)
    
    ##################################
    # 3D VARIABLES :
    ##################################
    
    print ('I am working on 3D vars..')
    
    
    # ===============================
    # Loop on 3D VARS :
    # ===============================
    
    for idx_3d in range(0,len(field_3d_name)):
        var_3d=field_3d_name[idx_3d]
        var_3d_udm=field_3d_units[idx_3d]
    
        # ====================================
        # Loop on VERTICAL LEVELS of 3D vars :
        # ====================================
        for idx_lev in range(0,len(field_3d_lev)):
            vlev=field_3d_lev[idx_lev]
            vlev_val=field_3d_lev_val[idx_lev]    

            # ===============================
            # Loop on PERIOD/SEASONS :
            # ===============================
            for idx_dt in range(0,len(period_of_analysis)):
                dt_lab=period_of_analysis[idx_dt]
    
                print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
                print ('# VERTICAL LEVEL: ',vlev)
                print ('# PERIOD: ',dt_lab)
    
                # Build the path/name of the nc file and open it 
                nc2open=model_path+model_fileprename+'3D_'+dt_lab+'_'+'allv'+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.nc'
                print ('Input file = ',nc2open)
                model = NC.Dataset(nc2open,'r')

                # Read lat, lon and fild values 
                if grid == 'uv2t':
                   # If destaggered with Max script:
                   lons = model.variables['lon'][:]
                   lats = model.variables['lat'][:]
                   # If destaggered with CDO:
                   #lons = model.variables['nav_lon'][:]
                   #lats = model.variables['nav_lat'][:]
                   #grid = 'T'
                   #T_threshold_3d=1
                   #thresh_sup=0.5
                   #thresh_inf=0.0
                   #print ('Grid: ',grid)
                   #
                   vals = np.sqrt(np.power(model.variables['ut'][:],2)+np.power(model.variables['vt'][:],2))
                   vals_du = model.variables['ut'][:]
                   vals_dv = model.variables['vt'][:]
                   # Sel vertical levs
                   print('vlev=',vlev_val)
                   vals=np.squeeze(vals[:,int(vlev_val),:,:])
                   vals_du=np.squeeze(vals_du[:,int(vlev_val),:,:])
                   vals_dv=np.squeeze(vals_dv[:,int(vlev_val),:,:])
                   vals_dv=np.ma.filled(vals_dv,0.000)
                   vals_du=np.ma.filled(vals_du,0.0000)
                else: 
                   lons = model.variables['nav_lon'][:]
                   lats = model.variables['nav_lat'][:]
                   vals = model.variables[var_3d][:]
                   # Sel vertical levs
                   vals=np.squeeze(vals[:,int(vlev_val),:,:])

                lats_red=[]
                lons_red=[]
                vals_red=[]
                # mask land and fix max and/or min values 
           #TMP     #thresh_inf = field_3d_inf[idx_3d]
           #TMP     #thresh_sup = field_3d_sup[idx_3d]
                thresh = 0.0000
                mask = np.abs(vals) == thresh 
                vals_ma = np.ma.masked_where(mask, vals)
                vals_ma_nan = np.ma.filled(vals_ma, -999.999)
                mask_rev= np.abs(vals) != thresh
                vals_ma_rev = np.ma.masked_where(mask_rev, vals)                
       
                # Plot the map and save in the path/name
                plotname=model_path+model_fileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.jpg'  
                print ('Plot path/name: ',plotname)
    
                plt.figure(figsize=(14,7)) #20,10
                plt.rc('font', size=12)
                plt.rcParams['figure.facecolor'] = 'black'
                # Plot the frame to the map
                plt.rcParams["axes.linewidth"]  = 1.5
                # Plot Title
                if dt_lab == 'yr':
                   plt.title ('Model '+model_label+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                else:
                   plt.title ('Model '+model_label+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')')
                # Read the coordinates for the plot 
                lon_0 = lons.mean()
                llcrnrlon = lons.min()
                urcrnrlon = lons.max()
                lat_0 = lats.mean()
                llcrnrlat = lats.min()
                urcrnrlat = lats.max()

                # Native vars case
                if grid != 'uv2t':
                   # Create the map
                   m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
                   xi, yi = m(lons, lats)
                   # Contour of land with black line
                   if grid == 'T':
                      land_value = [np.amin(vals_ma)]
                      contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
                   elif grid == 'U' or grid == 'V':
                      m.pcolor(xi,yi,np.squeeze(vals_ma_rev),cmap='Oranges')
                   # Plot the map
                   if grid == 'T':
                      # Fix theresholds on min and max values
                      if T_threshold_3d == 1:
                         cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=thresh_inf,vmax=thresh_sup)
                      else:
                         cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet')
                   elif grid == 'U' or grid == 'V':
                      cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=thresh_inf,vmax=thresh_sup)
                   # Add the grid
                   m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
                   m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
                   # Plot the legend and its label
                   cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
                   bar_label_string=var_3d+' ['+var_3d_udm+']'
                   cbar.set_label(bar_label_string)
                   vals_max=np.amax(vals_ma)
                   vals_min=np.amin(vals_ma)
                   text_min_x,text_min_y= m(-17,29.0)
                   text_max_x,text_max_y= m(32,29.0)
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)

                # Derived var (i) case
                else:
                   # Plot the map
                   plt.xlim(llcrnrlon,urcrnrlon)
                   plt.ylim(llcrnrlat,urcrnrlat)
                   plt.xlabel('LON')
                   plt.ylabel('LAT')
                   
                   vals_max=np.max(vals_ma)
                   thresh_sup=vals_max/2.0
                   # TO BE REMOVED (added for salinity):
                   if vlev == 10:
                      thresh_sup=0.225
                   elif vlev == 200:
                      thresh_sup=0.225

                   cs = plt.pcolor(lons,lats,np.squeeze(vals_ma),cmap='jet',vmin=0.0,vmax=thresh_sup) # vmin=0.0,vmax=thresh_sup
                   print ('Plot arrows (dir)..') #(1, 141, 380, 1307) 
                   for arr_idx_x in range(0,len(lons),10):
                       for arr_idx_y in range (0,len(lats),10):
                           if vals_du[arr_idx_y,arr_idx_x]!=0 and vals_dv[arr_idx_y,arr_idx_x]!=0:
                              arw = plt.arrow(lons[arr_idx_x],lats[arr_idx_y],vals_du[arr_idx_y,arr_idx_x]*3,vals_dv[arr_idx_y,arr_idx_x]*3,head_width=0.1 )
                   # Add the grid
                   plt.grid()
                   # Plot the legend and its label
                   land_value = -999.999
                   plt.colorbar(mappable=None,label='Currents'+' ['+var_3d_udm+']',orientation='horizontal',extend='max',aspect=60)
                   # Add land contour
                   if grid != 'uv2t':
                      contour = plt.contour(lons,lats,np.squeeze(vals_ma_nan),[land_value-500,land_value+500],colors='black')
                   else:
                      contour = plt.contour(lons,lats,np.squeeze(vals_ma_nan),[land_value-500,0.00],colors='black')
                   vals_max=np.amax(vals_ma)
                   vals_min=np.amin(vals_ma)
                   #text_min_x,text_min_y= m(-17,29.0)
                   #text_max_x,text_max_y= m(32,29.0)
                   text_max_x=32.0
                   text_max_y=29.0
                   text_min_x=-17.0
                   text_min_y=29.0
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)

                # Save and close 
                plt.savefig(plotname)
                plt.clf()
    
                # --------------------------------------------------
                # Plot the 3 map in different Med areas for currents
                if sub_plot_flag == 1:
                #if grid == 'U' or grid == 'V' or grid == 'T':

                   plotname=model_path+model_fileprename+'3D_sub_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.jpg'
                   print ('Plot path/name: ',plotname)

                   plt.figure(figsize=(20,12))
                   plt.rc('font', size=12)
                   plt.rcParams['figure.facecolor'] = 'black'
                   # Plot the frame to the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   # Plot Title
                   if dt_lab == 'yr':
                      plt.suptitle ('Model '+model_label+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                   else:
                      plt.suptitle ('Model '+model_label+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')')
       
                   # Loop on subareas
                   for idx_area in range (0,len(Area_name_ar)):
                       # Read the coordinates for the plot 
                       #lon_0 = lons.mean()
                       llcrnrlon = Area_minlon_ar[idx_area]
                       urcrnrlon = Area_maxlon_ar[idx_area]
                       #lat_0 = lats.mean()
                       llcrnrlat = Area_minlat_ar[idx_area]
                       urcrnrlat = Area_maxlat_ar[idx_area]
                       # Create the map
                       plt.subplot(2,2,idx_area+1)
                       plt.title (Area_name_ar[idx_area]+' Area')
                       m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc')
                       xi, yi = m(lons, lats)
                       # Plot the land and the map
                       m.pcolor(xi,yi,np.squeeze(vals_ma_rev),cmap='Oranges')
                       cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='gist_rainbow',vmin=Area_vmin_ar[idx_area],vmax=Area_vmax_ar[idx_area]) #jet
                       # Add the grid
                       m.drawparallels(np.arange(llcrnrlat,urcrnrlat, 1.), labels=[1,0,0,0], fontsize=10)
                       m.drawmeridians(np.arange(llcrnrlon,urcrnrlon, 2.), labels=[0,0,0,1], fontsize=10)
                       # Plot the legend and its label
                       cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
                       bar_label_string=var_3d+' ['+var_3d_udm+']'
                       cbar.set_label(bar_label_string)
                       vals_max=np.amax(vals_ma)
                       vals_min=np.amin(vals_ma)
                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()
    
    
    
    
    ##################################
    # 2D VARIABLES :
    ##################################
    
    print ('I am working on 2D vars..')
    
    
    # ===============================
    # Loop on 2D VARS :
    # ===============================
    
    for idx_2d in range(0,len(field_2d_name)):
        var_2d=field_2d_name[idx_2d]
        var_2d_udm=field_2d_units[idx_2d]
    
        # ===============================
        # Loop on PERIOD/SEASONS :
        # ===============================
        for idx_dt in range(0,len(period_of_analysis)):
            dt_lab=period_of_analysis[idx_dt]
    
            print ('# 2D VAR [udm]: ',var_2d,' [',var_2d_udm,']')
            print ('# PERIOD: ',dt_lab)
    
            # Build the path/name of the nc file and open it 
            if dataset_name != 'bathy':
               nc2open=model_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname+'.nc'
            else:
               nc2open=model_path+model_fileprename+'_'+model_postname+'.nc'

            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon and fild values 
            if grid == 'uv2t':
               lons = model.variables['lon'][:]
               lats = model.variables['lat'][:]
            else:
               lons = model.variables['nav_lon'][:]
               lats = model.variables['nav_lat'][:]
            vals = model.variables[var_2d][:]

            lats_red=[]
            lons_red=[]
            vals_red=[]
            # mask land and fix max and/or min values 
            if T_threshold_2d != 0:
               thresh_inf = field_2d_inf[idx_2d]
               thresh_sup = field_2d_sup[idx_2d]
            thresh = 0.0000
            #mask_sup = np.abs(vals) <= thresh_inf 
            #mask_inf = np.abs(vals) >= thresh_sup
            #vals_ma = np.ma.masked_where(mask_sup, vals)
            #vals_ma = np.ma.masked_where(mask_inf, vals_ma)
            mask = np.abs(vals) == thresh
            vals_ma = np.ma.masked_where(mask, vals)
    
            # Plot the map and save in the path/name
            if dataset_name != 'bathy':    
               plotname=workdir_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname+'.jpg'

               print ('Plot path/name: ',plotname)

               plt.figure(figsize=(20,10))
               plt.rc('font', size=12)
               # Plot Title
               if dt_lab == 'yr':
                  plt.title ('Model '+model_label+' MEAN '+var_2d+' ['+var_2d_udm+']'+' - DEPTH: 0 [m] - DT: '+dates_label)
               else:
                  plt.title ('Model '+model_label+' MEAN: '+var_2d+' ['+var_2d_udm+']'+' - DEPTH: 0 [m] - DT: '+dt_lab+' ('+dates_label+')')
               # Read the coordinates for the plot 
               lon_0 = lons.mean()
               llcrnrlon = lons.min()
               urcrnrlon = lons.max()
               lat_0 = lats.mean()
               llcrnrlat = lats.min()
               urcrnrlat = lats.max()
               print ('PROVA: ', llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)
               # Create the map
               m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
               xi, yi = m(lons, lats)
               # Plot the frame to the map
               plt.rcParams["axes.linewidth"]  = 1.25
               # Contour of land with black line
               #land_value = [0.0000]
               #if grid == 'T':
               #   contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
               m.drawmapboundary(fill_color='gray')
               # Plot the map
               if grid == 'T':
                  if T_threshold_2d != 0:
                     cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=thresh_inf) # ,vmin=thresh_inf
                  else:
                     vals_max=np.amax(vals_ma)
                     vals_min=np.amin(vals_ma)
                     thresh_inf=(abs(vals_max)+abs(vals_min))/2
                     cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=-0.25-0.20,vmax=-0.25+0.20)
               elif grid == 'U' or grid == 'V':
                    cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet')
               # Add the grid
               m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
               m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
               # Plot the legend and its label
               cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
               bar_label_string=var_2d+' ['+var_2d_udm+']'
               cbar.set_label(bar_label_string)
               vals_max=np.amax(vals_ma)
               vals_min=np.amin(vals_ma)
               text_min_x,text_min_y= m(-17,29.0)
               text_max_x,text_max_y= m(32,29.0)
               plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
               plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)
               # Save and close 
               plt.savefig(plotname)
               plt.clf()

            # Bathymetry case:
            else:
               if transpsect_flag == 1:
                  plotname=workdir_path+'transp_sections.jpg'
               elif tidegauge_flag == 1:
                  plotname=workdir_path+'tidegauges.jpg'
               else:
                  plotname=workdir_path+'model_bathymetry.jpg'

               print ('Plot path/name: ',plotname)
               plt.figure(figsize=(20,10)) # 20,10
               plt.rc('font', size=12) # 14
               # Plot Title
               if transpsect_flag == 1:
                  plt.title ('Bathymetry'+' ['+var_2d_udm+'] and transport sections')
               elif tidegauge_flag == 1:
                  plt.title ('Bathymetry'+' ['+var_2d_udm+'] and Tide-Gauges location')
               else:
                  plt.title ('Bathymetry'+' ['+var_2d_udm+']')
               # Read the coordinates for the plot 
               lon_0 = lons.mean()
               llcrnrlon = lons.min()
               llcrnrlon = -10.0
               urcrnrlon = lons.max()
               lat_0 = lats.mean()
               llcrnrlat = lats.min()
               urcrnrlat = lats.max()
               # Create the map
               m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
               xi, yi = m(lons, lats)
               # Plot the frame to the map
               plt.rcParams["axes.linewidth"]  = 1.25
               # Contour of land with black line
               land_value = [0.0000]
               if grid == 'T':
                  contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
               # Plot the map
               if grid == 'T':
                  if tidegauge_flag == 0:
                     cmap = mpl.cm.Blues(np.linspace(0,1,20))
                     cmap = mpl.colors.ListedColormap(cmap[5:,:-1])
                     cmap =  cmap.reversed()
                     cs = m.pcolor(xi,yi,-np.squeeze(vals_ma),cmap=cmap,vmax=thresh_inf,vmin=thresh_sup) # cmap='Blues_r',,vmin=thresh_inf,extend='min, sns.palplot(sns.dark_palette("blue"))
                  else:
                     cmap = mpl.cm.Greys(np.linspace(0,1,20))
                     cmap = mpl.colors.ListedColormap(cmap[:10,:-1])
                     cmap =  cmap.reversed()
                     cs = m.pcolor(xi,yi,-np.squeeze(vals_ma),cmap=cmap,vmax=thresh_inf,vmin=thresh_sup) 
               elif grid == 'U' or grid == 'V':
                    cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet')
               # Add the grid
               m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
               m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
               # Plot the legend and its label
               cbar = m.colorbar(cs, location='right', pad="10%")
               bar_label_string=var_2d+' ['+var_2d_udm+']'
               cbar.set_label(bar_label_string)
               if tidegauge_flag == 1:
                 # ADD TG LOCATIONS
                 xp, yp = m(13.50833333333333333333,43.62694444444444444443)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#0000ff',label='28') # TG ancona 
                 plt.text(xp,yp,'28', fontsize=12,backgroundcolor='#0000ff',alpha=1,color='white')
                 xp, yp = m(8.30500000000000000000,39.13555555555555555555)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='19') # TG carloforte 
                 plt.text(xp,yp,'19', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(15.08472222222222222221,37.49138888888888888888)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ffa500',label='32') # TG catania 
                 plt.text(xp,yp,'32', fontsize=12,backgroundcolor='#ffa500',alpha=1,color='black')
                 xp, yp = m(8.01694444444444444443,43.88111111111111111110)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='18') # TG imperia 
                 plt.text(xp,yp,'18', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(12.61000000000000000000,35.49138888888888888888)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='27') # TG lampedusa 
                 plt.text(xp,yp,'27', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(10.28805555555555555555,43.54222222222222222221)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='25') # TG livorno 
                 plt.text(xp,yp,'25', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(15.55916666666666666666,38.18638888888888888888)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ffa500',label='35') # TG messina 
                 plt.text(xp,yp,'35,36', fontsize=12,backgroundcolor='#ffa500',alpha=1,color='black')
                 xp, yp = m(14.40666666666666666666,42.35583333333333333333)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#0000ff',label='31') # TG ortona 
                 plt.text(xp,yp,'31', fontsize=12,backgroundcolor='#0000ff',alpha=1,color='white')
                 xp, yp = m(15.27111111111111111110,40.01694444444444444443)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='34') # TG palinuro 
                 plt.text(xp,yp,'34', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(13.52527777777777777777,37.28805555555555555555)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='29') # TG portoempedocle 
                 plt.text(xp,yp,'29', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 ##xp, yp = m(8.76290,41.92270)
                 ###plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='21') # TG AjaccioTG 
                 ##plt.text(xp,yp,'21', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(8.40666666666666666666,40.84722222222222222221)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='20') # TG portotorres 
                 plt.text(xp,yp,'20', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 #xp, yp = m(15.64388888888888888888,38.11861111111111111110)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ffa500',label='36') # TG reggiocalabria 
                 #plt.text(xp,yp,'36', fontsize=12,backgroundcolor='#ffa500',alpha=1,color='black')
                 xp, yp = m(13.76250000000000000000,45.64388888888888888888)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#0000ff',label='30') # TG trieste 
                 plt.text(xp,yp,'30', fontsize=12,backgroundcolor='#0000ff',alpha=1,color='white')
                 xp, yp = m(12.42361111111111111110,45.42361111111111111110)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#0000ff',label='26') # TG venezia 
                 plt.text(xp,yp,'26', fontsize=12,backgroundcolor='#0000ff',alpha=1,color='white')
                 xp, yp = m(16.16944444444444444443,41.89805555555555555555)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#0000ff',label='37') # TG vieste 
                 plt.text(xp,yp,'37', fontsize=12,backgroundcolor='#0000ff',alpha=1,color='white')
                 xp, yp = m(8.66290,41.92270) # 8.76290
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='21') # TG AjaccioTG 
                 plt.text(xp,yp,'21', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 #xp, yp = m(-5.39833,36.1769)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff0000',label='2') # TG AlgecirasTG 
                 #plt.text(xp,yp,'2', fontsize=12,backgroundcolor='#ff0000',alpha=1,color='black')
                 xp, yp = m(-2.478,36.83)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='6') # TG AlmeriaTG 
                 plt.text(xp,yp,'6', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(2.163,41.342)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='10') # TG BarcelonaTG 
                 plt.text(xp,yp,'10', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(9.34983,42.96578)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='23') # TG CenturiTG 
                 plt.text(xp,yp,'23', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(1.44972,38.9111)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='9') # TG IbizaTG 
                 plt.text(xp,yp,'9', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(8.93524,42.63960)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='22') # TG IleRousseTG 
                 plt.text(xp,yp,'22', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 ##xp, yp = m(6.93377,43.48353)
                 ###plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='16') # TG LaFigueiretteTG 
                 ##plt.text(xp,yp,'16', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(-4.417,36.712)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff0000',label='3') # TG MalagaTG 
                 plt.text(xp,yp,'3', fontsize=12,backgroundcolor='#ff0000',alpha=1,color='black')
                 xp, yp = m(5.35370,43.27850)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='15') # TG MarseilleTG 
                 plt.text(xp,yp,'15', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(-2.918,35.291)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff0000',label='5') # TG MelillaTG 
                 plt.text(xp,yp,'5', fontsize=12,backgroundcolor='#ff0000',alpha=1,color='black')
                 xp, yp = m(7.42370,43.73300)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='17') # TG MonacoTG 
                 plt.text(xp,yp,'17', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(6.93377,43.48353)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='16') # TG LaFigueiretteTG 
                 plt.text(xp,yp,'16', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(-3.524,36.72)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff0000',label='4') # TG MotrilTG 
                 plt.text(xp,yp,'4', fontsize=12,backgroundcolor='#ff0000',alpha=1,color='black')
                 ##xp, yp = m(3.06410,43.01471)
                 ###plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='12') # TG PortLaNouvelleTG 
                 ##plt.text(xp,yp,'12', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(2.6375,39.5603)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='11') # TG PalmadeMallorcaTG 
                 plt.text(xp,yp,'11', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 ##xp, yp = m(3.06410,43.01471)
                 ###plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='12') # TG PortLaNouvelleTG 
                 ##plt.text(xp,yp,'12', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(3.10730,42.52010)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='13') # TG PortVendresTG 
                 plt.text(xp,yp,'13', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 #xp, yp = m(-0.206,39.634)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='8') # TG SaguntoTG 
                 #plt.text(xp,yp,'8', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(3.70170,43.40000)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='14') # TG SeteTG 
                 plt.text(xp,yp,'14', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(3.06410,43.01471)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='12') # TG PortLaNouvelleTG 
                 plt.text(xp,yp,'12', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(9.40383,41.85686)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='24') # TG SolenzaraTG 
                 plt.text(xp,yp,'24', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(-5.60361,36.0064)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff0000',label='1') # TG TarifaTG 
                 plt.text(xp,yp,'1,2', fontsize=12,backgroundcolor='#ff0000',alpha=1,color='black')
                 xp, yp = m(-0.33,39.46)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='7') # TG ValenciaTG 
                 plt.text(xp,yp,'7,8', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')
                 xp, yp = m(33.340228,34.726315)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff00ff',label='38') # TG zygi1TG 
                 plt.text(xp,yp,'38', fontsize=12,backgroundcolor='#ff00ff',alpha=1,color='black')
                 xp, yp = m(36.176769256592,36.594230651855)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#ff00ff',label='39') # TG iske 
                 plt.text(xp,yp,'39', fontsize=12,backgroundcolor='#ff00ff',alpha=1,color='black')
                 xp, yp = m(15.1933,38.784)
                 #plt.scatter(xp, yp, s=100, alpha=1, c='#32cd32',label='33') # TG Ginostra 
                 plt.text(xp,yp,'33', fontsize=12,backgroundcolor='#32cd32',alpha=1,color='black')


               if transpsect_flag == 1:

#                       #TRA_Gibraltar_6_15
#                       lonsss1,latsss1=m(-6.15,35.0)
#                       lonsss2,latsss2=m(-6.15,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#
#                       #TRA_Gibraltar_6_10
#                       lonsss1,latsss1=m(-6.10,35.0)
#                       lonsss2,latsss2=m(-6.10,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_6_06
#                       lonsss1,latsss1=m(-6.06,35.0)
#                       lonsss2,latsss2=m(-6.06,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_6_02
#                       lonsss1,latsss1=m(-6.02,35.0)
#                       lonsss2,latsss2=m(-6.02,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_6_98
#                       lonsss1,latsss1=m(-5.98,35.0)
#                       lonsss2,latsss2=m(-5.98,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Gibraltar_5_94
                       lonsss1,latsss1=m(-5.94,35.5)
                       lonsss2,latsss2=m(-5.94,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
 
#                       #TRA_Gibraltar_5_90
#                       lonsss1,latsss1=m(-5.90,35.0)
#                       lonsss2,latsss2=m(-5.90,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_85
#                       lonsss1,latsss1=m(-5.85,35.0)
#                       lonsss2,latsss2=m(-5.85,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_81
#                       lonsss1,latsss1=m(-5.81,35.0)
#                       lonsss2,latsss2=m(-5.81,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_77
#                       lonsss1,latsss1=m(-5.77,35.5)
#                       lonsss2,latsss2=m(-5.77,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_73
#                       lonsss1,latsss1=m(-5.73,35.0)
#                       lonsss2,latsss2=m(-5.73,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_69
#                       lonsss1,latsss1=m(-5.69,35.0)
#                       lonsss2,latsss2=m(-5.69,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_65
#                       lonsss1,latsss1=m(-5.65,35.0)
#                       lonsss2,latsss2=m(-5.65,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_60
#                       lonsss1,latsss1=m(-5.60,35.0)
#                       lonsss2,latsss2=m(-5.60,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_56
#                       lonsss1,latsss1=m(-5.56,35.0)
#                       lonsss2,latsss2=m(-5.56,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_52
#                       lonsss1,latsss1=m(-5.52,35.0)
#                       lonsss2,latsss2=m(-5.52,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Gibraltar_5_48
                       lonsss1,latsss1=m(-5.48,35.5)
                       lonsss2,latsss2=m(-5.48,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_44
#                       lonsss1,latsss1=m(-5.44,35.0)
#                       lonsss2,latsss2=m(-5.44,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_40
#                       lonsss1,latsss1=m(-5.40,35.0)
#                       lonsss2,latsss2=m(-5.40,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_35
#                       lonsss1,latsss1=m(-5.35,35.0)
#                       lonsss2,latsss2=m(-5.35,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Gibraltar_5_31
                       lonsss1,latsss1=m(-5.31,35.5)
                       lonsss2,latsss2=m(-5.31,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_27
#                       lonsss1,latsss1=m(-5.27,35.0)
#                       lonsss2,latsss2=m(-5.27,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_23
#                       lonsss1,latsss1=m(-5.23,35.0)
#                       lonsss2,latsss2=m(-5.23,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_19
#                       lonsss1,latsss1=m(-5.19,35.0)
#                       lonsss2,latsss2=m(-5.19,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_15
#                       lonsss1,latsss1=m(-5.15,35.0)
#                       lonsss2,latsss2=m(-5.15,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_10
#                       lonsss1,latsss1=m(-5.10,35.0)
#                       lonsss2,latsss2=m(-5.10,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_06
#                       lonsss1,latsss1=m(-5.06,35.0)
#                       lonsss2,latsss2=m(-5.06,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Gibraltar_5_02
#                       lonsss1,latsss1=m(-5.02,35.0)
#                       lonsss2,latsss2=m(-5.02,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Sicily
                       lonsss1,latsss1=m(9.91,36.95)
                       lonsss2,latsss2=m(14.46,36.96)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
 
                       #TRA_Otranto
                       lonsss1,latsss1=m(18.04,40.25)
                       lonsss2,latsss2=m(19.7,40.26)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
 
                       #TRA_Corsica
                       lonsss1,latsss1=m(9.12,42.5)
                       lonsss2,latsss2=m(11.58,42.51)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_OBL
#                       lonsss1,latsss1=m(15.54,38.2)
#                       lonsss2,latsss2=m(15.75,38.21)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_00
#                       lonsss1,latsss1=m(15.30,38.0)
#                       lonsss2,latsss2=m(15.90,38.00)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_04
#                       lonsss1,latsss1=m(15.30,38.04)
#                       lonsss2,latsss2=m(15.90,38.04)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_08
#                       lonsss1,latsss1=m(15.30,38.08)
#                       lonsss2,latsss2=m(15.90,38.08)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_13
#                       lonsss1,latsss1=m(15.30,38.13)
#                       lonsss2,latsss2=m(15.90,38.13)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_17
#                       lonsss1,latsss1=m(15.30,38.17)
#                       lonsss2,latsss2=m(15.90,38.17)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Messina_21
                       lonsss1,latsss1=m(15.55,38.21)
                       lonsss2,latsss2=m(15.90,38.21)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Messina_25
#                       lonsss1,latsss1=m(15.55,38.25)
#                       lonsss2,latsss2=m(15.90,38.25)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
# 
                       #TRA_Messina_29
                       lonsss1,latsss1=m(15.55,38.29)
                       lonsss2,latsss2=m(15.90,38.29)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
 
                       #TRA_Dardanelles
                       lonsss1,latsss1=m(26.16,39.97)
                       lonsss2,latsss2=m(26.16,40.11)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
# 
#                       #TRA_Ciprus
#                       lonsss1,latsss1=m(33.50,35.00)
#                       lonsss2,latsss2=m(36.00,35.00)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')

            # Save and close 
            plt.savefig(plotname)
            plt.clf()

            # --------------------------------------------------
            # Plot the bathy map in different Med areas
            if sub_plot_flag == 1 and grid != 'uv2t' and transpsect_flag == 1:

                   plotname=workdir_path+'transp_sections_sub.jpg'
                   print ('Plot path/name: ',plotname)

                   plt.figure(figsize=(20,12))
                   plt.rc('font', size=12)
                   plt.rcParams['figure.facecolor'] = 'black'
                   # Plot the frame to the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   # Plot Title
                   if dt_lab == 'yr':
                      plt.suptitle ('Bathymetry'+' ['+var_2d_udm+'] and Transport Sections -- Med Sub-regions')
                   # Loop on subareas
                   for idx_area in range (0,len(Area_name_ar)):
                       # Read the coordinates for the plot 
                       #lon_0 = lons.mean()
                       llcrnrlon = Area_minlon_ar[idx_area]
                       urcrnrlon = Area_maxlon_ar[idx_area]
                       #lat_0 = lats.mean()
                       llcrnrlat = Area_minlat_ar[idx_area]
                       urcrnrlat = Area_maxlat_ar[idx_area]
                       # Create the map
                       plt.subplot(2,2,idx_area+1)
                       plt.title (Area_name_ar[idx_area]+' Area')
                       m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc')
                       xi, yi = m(lons, lats)
                       # Plot the land and the map
                       cmap = mpl.cm.Blues(np.linspace(0,1,20))
                       cmap = mpl.colors.ListedColormap(cmap[5:,:-1])
                       cmap =  cmap.reversed()
                       cs = m.pcolor(xi,yi,-np.squeeze(vals_ma),cmap=cmap,vmax=Area_vmin_ar[idx_area],vmin=Area_vmax_ar[idx_area]) #cmap='gist_rainbow'
                       # Contour of land with black line
                       land_value = [-0]
                       if grid == 'T':
                          contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
                       # Add the grid
                       m.drawparallels(np.arange(llcrnrlat,urcrnrlat, 1.), labels=[1,0,0,0], fontsize=10)
                       m.drawmeridians(np.arange(llcrnrlon,urcrnrlon, 2.), labels=[0,0,0,1], fontsize=10)
                       # Plot the legend and its label
                       cbar = m.colorbar(cs, location='right', pad="10%")
                       bar_label_string=var_2d+' ['+var_2d_udm+']'
                       cbar.set_label(bar_label_string)

#                       #TRA_Gibraltar_6_15
#                       lonsss1,latsss1=m(-6.15,35.0)
#                       lonsss2,latsss2=m(-6.15,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_6_10
#                       lonsss1,latsss1=m(-6.10,35.0)
#                       lonsss2,latsss2=m(-6.10,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_6_06
#                       lonsss1,latsss1=m(-6.06,35.0)
#                       lonsss2,latsss2=m(-6.06,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_6_02
#                       lonsss1,latsss1=m(-6.02,35.0)
#                       lonsss2,latsss2=m(-6.02,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_6_98
#                       lonsss1,latsss1=m(-5.98,35.0)
#                       lonsss2,latsss2=m(-5.98,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Gibraltar_5_94
                       lonsss1,latsss1=m(-5.94,35.5)
                       lonsss2,latsss2=m(-5.94,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
                       
#                       #TRA_Gibraltar_5_90
#                       lonsss1,latsss1=m(-5.90,35.0)
#                       lonsss2,latsss2=m(-5.90,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_85
#                       lonsss1,latsss1=m(-5.85,35.0)
#                       lonsss2,latsss2=m(-5.85,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_81
#                       lonsss1,latsss1=m(-5.81,35.0)
#                       lonsss2,latsss2=m(-5.81,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_77
#                       lonsss1,latsss1=m(-5.77,35.5)
#                       lonsss2,latsss2=m(-5.77,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_73
#                       lonsss1,latsss1=m(-5.73,35.0)
#                       lonsss2,latsss2=m(-5.73,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_69
#                       lonsss1,latsss1=m(-5.69,35.0)
#                       lonsss2,latsss2=m(-5.69,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_65
#                       lonsss1,latsss1=m(-5.65,35.0)
#                       lonsss2,latsss2=m(-5.65,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_60
#                       lonsss1,latsss1=m(-5.60,35.0)
#                       lonsss2,latsss2=m(-5.60,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_56
#                       lonsss1,latsss1=m(-5.56,35.0)
#                       lonsss2,latsss2=m(-5.56,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_52
#                       lonsss1,latsss1=m(-5.52,35.0)
#                       lonsss2,latsss2=m(-5.52,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Gibraltar_5_48
                       lonsss1,latsss1=m(-5.48,35.5)
                       lonsss2,latsss2=m(-5.48,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_44
#                       lonsss1,latsss1=m(-5.44,35.0)
#                       lonsss2,latsss2=m(-5.44,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_40
#                       lonsss1,latsss1=m(-5.40,35.0)
#                       lonsss2,latsss2=m(-5.40,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_35
#                       lonsss1,latsss1=m(-5.35,35.0)
#                       lonsss2,latsss2=m(-5.35,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Gibraltar_5_31
                       lonsss1,latsss1=m(-5.31,35.5)
                       lonsss2,latsss2=m(-5.31,36.5)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
                       
#                       #TRA_Gibraltar_5_27
#                       lonsss1,latsss1=m(-5.27,35.0)
#                       lonsss2,latsss2=m(-5.27,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_23
#                       lonsss1,latsss1=m(-5.23,35.0)
#                       lonsss2,latsss2=m(-5.23,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_19
#                       lonsss1,latsss1=m(-5.19,35.0)
#                       lonsss2,latsss2=m(-5.19,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_15
#                       lonsss1,latsss1=m(-5.15,35.0)
#                       lonsss2,latsss2=m(-5.15,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_10
#                       lonsss1,latsss1=m(-5.10,35.0)
#                       lonsss2,latsss2=m(-5.10,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_06
#                       lonsss1,latsss1=m(-5.06,35.0)
#                       lonsss2,latsss2=m(-5.06,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Gibraltar_5_02
#                       lonsss1,latsss1=m(-5.02,35.0)
#                       lonsss2,latsss2=m(-5.02,36.5)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Sicily
                       lonsss1,latsss1=m(9.91,36.95)
                       lonsss2,latsss2=m(14.46,36.96)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
                       
                       #TRA_Otranto
                       lonsss1,latsss1=m(18.04,40.25)
                       lonsss2,latsss2=m(19.7,40.26)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
                       
                       #TRA_Corsica
                       lonsss1,latsss1=m(9.12,42.5)
                       lonsss2,latsss2=m(11.58,42.51)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_OBL
#                       lonsss1,latsss1=m(15.54,38.2)
#                       lonsss2,latsss2=m(15.75,38.21)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_00
#                       lonsss1,latsss1=m(15.30,38.0)
#                       lonsss2,latsss2=m(15.90,38.00)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_04
#                       lonsss1,latsss1=m(15.30,38.04)
#                       lonsss2,latsss2=m(15.90,38.04)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_08
#                       lonsss1,latsss1=m(15.30,38.08)
#                       lonsss2,latsss2=m(15.90,38.08)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_13
#                       lonsss1,latsss1=m(15.30,38.13)
#                       lonsss2,latsss2=m(15.90,38.13)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_17
#                       lonsss1,latsss1=m(15.30,38.17)
#                       lonsss2,latsss2=m(15.90,38.17)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Messina_21
                       lonsss1,latsss1=m(15.55,38.21)
                       lonsss2,latsss2=m(15.90,38.21)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Messina_25
#                       lonsss1,latsss1=m(15.55,38.25)
#                       lonsss2,latsss2=m(15.90,38.25)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
#                       
                       #TRA_Messina_29
                       lonsss1,latsss1=m(15.55,38.29)
                       lonsss2,latsss2=m(15.90,38.29)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
                       
                       #TRA_Dardanelles
                       lonsss1,latsss1=m(26.16,39.97)
                       lonsss2,latsss2=m(26.16,40.11)
                       lonsss=[lonsss1,lonsss2]
                       latsss=[latsss1,latsss2]
                       plt.plot(lonsss,latsss,color='red')
#                       
#                       #TRA_Ciprus
#                       lonsss1,latsss1=m(33.50,35.00)
#                       lonsss2,latsss2=m(36.00,35.00)
#                       lonsss=[lonsss1,lonsss2]
#                       latsss=[latsss1,latsss2]
#                       plt.plot(lonsss,latsss,color='red')
         
            # Save and close 
            plt.savefig(plotname)
            plt.clf()


#########################################
######### DIFF BETWEEN MODELS ###########

elif num_of_models == 2:
    # ======================================
    # Diff between 1 and 2  model datasets :
    # ======================================
    model_dt_1=model_label[0]
    model_dt_2=model_label[1]
 
    ##################################
    # 3D VARIABLES :
    ##################################

    print ('I am working on 3D vars..')
    print ('Diff case...')

    # ===============================
    # Loop on 3D VARS :
    # ===============================

    for idx_3d in range(0,len(field_3d_name)):
        var_3d=field_3d_name[idx_3d]
        var_3d_udm=field_3d_units[idx_3d]

        # ====================================
        # Loop on VERTICAL LEVELS of 3D vars :
        # ====================================
        for idx_lev in range(0,len(field_3d_lev)):
            vlev=field_3d_lev[idx_lev]
            vlev_val=field_3d_lev_val[idx_lev]

            # ===============================
            # Loop on PERIOD/SEASONS :
            # ===============================
            for idx_dt in range(0,len(period_of_analysis)):
                dt_lab=period_of_analysis[idx_dt]

                print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
                print ('# VERTICAL LEVEL: ',vlev)
                print ('# PERIOD: ',dt_lab)

                # Build the path/name of the nc file and open it 
                nc2open1=model_path+model_infileprename+'3D_'+dt_lab+'_'+'allv'+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'.nc'
                nc2open2=model_path+model_infileprename+'3D_'+dt_lab+'_'+'allv'+'_'+var_3d+'_'+dates_label+'_'+model_postname[1]+'.nc'
                print ('Input files = ',nc2open1,nc2open2)
                model1 = NC.Dataset(nc2open1,'r')
                model2 = NC.Dataset(nc2open2,'r')

                # Read lat, lon and fild values 
                if grid == 'uv2t':
                   lons = model1.variables['lon'][:]
                   lats = model2.variables['lat'][:]
                   vals = model1.variables['i'][:]-model2.variables['i'][:]
                   #vals = np.sqrt(np.power(model1.variables['ut'][:]-model2.variables['ut'][:],2)+np.power(model1.variables['vt'][:]-model2.variables['vt'][:],2))
                   vals_du = model1.variables['ut'][:]-model2.variables['ut'][:]
                   vals_dv = model1.variables['vt'][:]-model2.variables['vt'][:]
                   # Sel vertical levs
                   print('vlev=',vlev_val)
                   vals=np.squeeze(vals[:,int(vlev_val),:,:])
                   vals_du=np.squeeze(vals_du[:,int(vlev_val),:,:])
                   vals_dv=np.squeeze(vals_dv[:,int(vlev_val),:,:])
                   vals_dv=np.ma.filled(vals_dv,0.000)
                   vals_du=np.ma.filled(vals_du,0.0000)
                   #print ('PROVA: ',vals.shape,vals_du.shape,vals_dv.shape)
                else:
                   lons = model1.variables['nav_lon'][:]
                   lats = model2.variables['nav_lat'][:]
                   vals = model1.variables[var_3d][:]-model2.variables[var_3d][:]
                   # Sel vertical levs
                   print('vlev=',vlev_val)
                   vals=np.squeeze(vals[:,int(vlev_val),:,:])

                lats_red=[]
                lons_red=[]
                vals_red=[]
                # Built the palette 
                vals_max=np.max(vals)
                vals_min=np.min(vals)

                # mask land and fix max and/or min values 
                thresh_inf = field_3d_inf[idx_3d]
                thresh_sup = field_3d_sup[idx_3d]
                thresh = 0.0000
                mask = np.abs(vals) == thresh
                vals_ma = np.ma.masked_where(mask, vals)
                vals_ma_nan = np.ma.filled(vals_ma, -999.999)
                mask_rev= np.abs(vals) != thresh
                vals_ma_rev = np.ma.masked_where(mask_rev, vals)

                # Built the palette 
                vals_max=np.amax(vals_ma)
                vals_min=np.amin(vals_ma)
                # auto defined threshold
                if Tdiff_threshold_3d == 0:
                   #vals_thresc=min(abs(vals_max),abs(vals_min))
                   if var_3d == 'vosaline':
                        vals_thresc=0.1
                   #   vals_thresc=min(abs(vals_max),abs(vals_min))
                   #   if str(vlev) == '1':
                   #      vals_thresc=0.5
                   #   if str(vlev) == '30':
                   #      vals_thresc=0.5 # 0.5
                   #   if str(vlev) == '150':
                   #      vals_thresc=0.5
                        #vals_thresc=min(abs(vals_max),abs(vals_min))-1.0
                        #vals_thresc=min(abs(vals_max),abs(vals_min))
                   else:
                         vals_thresc=0.6
                      #vals_thresc=(abs(vals_max)+abs(vals_min))/2
                   #   vals_thresc=min(abs(vals_max),abs(vals_min))+0.5                      
                   #   if str(vlev) == '1': 
                   #      vals_thresc=0.80 # 1.2
                   #   if str(vlev) == '30':
                   #      vals_thresc=0.80 #2.2
                   #   if str(vlev) == '150':
                   #      vals_thresc=0.80

                # user deined thresholds:
                elif Tdiff_threshold_3d == 1:
                     vals_thresc=tdiff_th
                print ('vals_thresc = [-',vals_thresc,',',vals_thresc,']')


                # Plot the map and save in the path/name
                plotname=model_path+model_fileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
                print ('Plot path/name: ',plotname)

                plt.figure(figsize=(14,7)) #20,10
                plt.rc('font', size=12)
                # Plot Title
                if dt_lab == 'yr':
                   plt.title ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                else:
                   plt.title ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')')
                # Read the coordinates for the plot 
                lon_0 = lons.mean()
                llcrnrlon = lons.min()
                urcrnrlon = lons.max()
                lat_0 = lats.mean()
                llcrnrlat = lats.min()
                urcrnrlat = lats.max()
                #
                # Native vars case
                if grid != 'uv2t':

                   # Create the map
                   m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
                   xi, yi = m(lons, lats)
                   # Plot the frame to the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   # Contour of land with black line
                   land_value = [0.0000]
                   # Plot the map
                   if var_3d == 'votemper':
                      cmap='bwr'
                   elif var_3d == 'vosaline':
                      cmap='BrBG'
                   else:
                      cmap='bwr'
                   cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap=cmap,vmin=-vals_thresc,vmax=vals_thresc)
                   text_min_x,text_min_y= m(-17,29.0)
                   text_max_x,text_max_y= m(32,29.0)
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)
                   # Add the grid
                   m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
                   m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
                   # Plot the legend and its label
                   cbar = m.colorbar(cs,location='bottom', pad="10%",extend='both') # location='bottom', pad="10%",extend='both'
                   bar_label_string=var_3d+' ['+var_3d_udm+']'
                   cbar.set_label(bar_label_string)
                   # Land contour
                   land_value = [-999.999,vals_min-100] #vals_min-100]
                   #m.contourf(xi,yi,np.squeeze(vals_ma_nan),land_value,colors='black')
                   m.drawmapboundary(fill_color='gray') # fill_color='gray'
                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()

                # Derived var (i-uv2t) case
                else:
                   #
                   # Plot the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   plt.xlim(llcrnrlon,urcrnrlon)
                   plt.ylim(llcrnrlat,urcrnrlat)
                   plt.xlabel('LON')
                   plt.ylabel('LAT')
                   ##### TO BE RM: 
                   #if vlev == 10:
                   #   vals_thresc=0.03 # 0.06 0.22; 0.15
                   #elif vlev == 200:
                   #   vals_thresc=0.03 #0.066 vale
                   #
                   cs = plt.pcolor(lons,lats,np.squeeze(vals_ma),cmap='RdBu_r',vmin=-vals_thresc,vmax=vals_thresc) # cmap='jet',vmin=0,vmax=vals_thresc
                   ###cs = plt.contourf(lons,lats,np.squeeze(vals_ma),[-0.042,-0.036,-0.030,-0.024,-0.018,-0.012,-0.006,0.00,0.006,0.012,0.018,0.024,0.030,0.036,0.042],cmap='seismic',extend='both')
                   # Comparison with Agresti 10m
                   #cs = plt.contourf(lons,lats,np.squeeze(vals_ma),[-0.02,0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20],cmap='jet')
                   ## Comparison with Agresti 10m
                   ##cs = plt.contourf(lons,lats,np.squeeze(vals_ma),[-0.006,0.00,0.006,0.012,0.018,0.024,0.030,0.036,0.042,0.048,0.054,0.060],cmap='jet',extend='both')
                   # 
                   # Add the grid
                   plt.grid()
                   # Plot the legend and its label
                   land_value = -999.999
                   plt.colorbar(mappable=None,label='Currents'+' ['+var_3d_udm+']',orientation='horizontal',extend='both',aspect=60)
                   # Add land contour
                   contour = plt.contourf(lons,lats,np.squeeze(vals_ma_nan),[-999.999,vals_min-100],colors='black')
                   print ('Plot arrows (dir)..') #(1, 141, 380, 1307) 
                   for arr_idx_x in range(0,len(lons),10): 
                       for arr_idx_y in range (0,len(lats),10):
                           if vals_du[arr_idx_y,arr_idx_x]!=0 and vals_dv[arr_idx_y,arr_idx_x]!=0:
                              arw = plt.arrow(lons[arr_idx_x],lats[arr_idx_y],vals_du[arr_idx_y,arr_idx_x]*3,vals_dv[arr_idx_y,arr_idx_x]*3,head_width=0.1 )
                   #,head_width=10,head_lenght=10)
                   #
                   vals_max=np.amax(vals_ma)
                   vals_min=np.amin(vals_ma)
                   #text_min_x,text_min_y= m(-17,29.0)
                   #text_max_x,text_max_y= m(32,29.0)
                   text_max_x=32.0
                   text_max_y=29.0
                   text_min_x=-17.0
                   text_min_y=29.0
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)

                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()


                # --------------------------------------------------
                # Plot the 3 map in different Med areas for currents
                if sub_plot_flag == 1 and grid != 'uv2t':
                #if grid == 'U' or grid == 'V' or grid == 'T':
               
                   plotname=model_path+model_fileprename+'3D_sub_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
                   print ('Plot path/name: ',plotname)

                   plt.figure(figsize=(20,12))
                   plt.rc('font', size=12)
                   plt.rcParams['figure.facecolor'] = 'black'
                   # Plot the frame to the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   # Plot Title
                   if dt_lab == 'yr':
                      plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                   else:
                      plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')')
                   # Loop on subareas
                   for idx_area in range (0,len(Area_name_ar)):
                       # Read the coordinates for the plot 
                       #lon_0 = lons.mean()
                       llcrnrlon = Area_minlon_ar[idx_area]
                       urcrnrlon = Area_maxlon_ar[idx_area]
                       #lat_0 = lats.mean()
                       llcrnrlat = Area_minlat_ar[idx_area]
                       urcrnrlat = Area_maxlat_ar[idx_area]
                       # Create the map
                       plt.subplot(2,2,idx_area+1)
                       plt.title (Area_name_ar[idx_area]+' Area')
                       m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc')
                       xi, yi = m(lons, lats)
                       # Plot the land and the map
                       m.pcolor(xi,yi,np.squeeze(vals_ma_rev),cmap='gray') #cmap='Oranges'
                       cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='bwr',vmin=Area_vmin_ar[idx_area],vmax=Area_vmax_ar[idx_area]) #cmap='gist_rainbow'
                       # Add the grid
                       m.drawparallels(np.arange(llcrnrlat,urcrnrlat, 1.), labels=[1,0,0,0], fontsize=10)
                       m.drawmeridians(np.arange(llcrnrlon,urcrnrlon, 2.), labels=[0,0,0,1], fontsize=10)
                       # Plot the legend and its label
                       cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
                       bar_label_string=var_3d+' ['+var_3d_udm+']'
                       cbar.set_label(bar_label_string)
                       vals_max=np.amax(vals_ma)
                       vals_min=np.amin(vals_ma)

                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()

                if sub_plot_flag == 1 and grid == 'uv2t':

                   plotname=model_path+model_fileprename+'3D_sub_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
                   print ('Plot path/name: ',plotname)

                   plt.figure(figsize=(10,10))
                   plt.rc('font', size=12)
                   plt.rcParams['figure.facecolor'] = 'black'
                   # Plot the frame to the map
                   plt.rcParams["axes.linewidth"]  = 1.25
                   # Plot Title
                   if dt_lab == 'yr':
                      plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                   else:
                      plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')')
                   # Loop on subareas
                   for idx_area in range (0,len(Area_name_ar)):
                   #idx_area=0
                       # Read the coordinates for the plot 
                       #lon_0 = lons.mean()
                       llcrnrlon = Area_minlon_ar[idx_area]
                       urcrnrlon = Area_maxlon_ar[idx_area]
                       #lat_0 = lats.mean()
                       llcrnrlat = Area_minlat_ar[idx_area]
                       urcrnrlat = Area_maxlat_ar[idx_area]
                       #
                   ## TO BE REMOVED
                   #llcrnrlon=-8
                   #urcrnrlon=-2
                   #llcrnrlat=34
                   #urcrnrlat=38

                   # Create the map
                       plt.subplot(2,2,idx_area+1)
                       plt.title (Area_name_ar[idx_area]+' Area')
                       if dt_lab == 'yr':
                          plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dates_label)
                       else:
                          plt.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - DEPTH: '+str(vlev)+'[m] - DT: '+dt_lab+' ('+dates_label+')') 
                   #
                       plt.xlim(llcrnrlon,urcrnrlon)
                       plt.ylim(llcrnrlat,urcrnrlat)
                       plt.xlabel('LON')
                       plt.ylabel('LAT')
                       #m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc')
                       #xi, yi = m(lons, lats)
                       # Plot the land and the map
                       #m.pcolor(lons,lats,np.squeeze(vals_ma_rev),cmap='Oranges')
                       cs = plt.pcolor(lons,lats,np.squeeze(vals_ma),cmap='RdBu_r',vmin=Area_vmin_ar[idx_area],vmax=Area_vmax_ar[idx_area])
                       # Plot the legend and its label
                       plt.colorbar(mappable=None,label='Currents'+' ['+var_3d_udm+']',orientation='horizontal',extend='both') #,aspect=60)
                       # Add the grid
                       plt.grid()
                       land_value = -999.999
                       # Add land contour
                       contour = plt.contourf(lons,lats,np.squeeze(vals_ma_nan),[-999.999,vals_min-100],colors='black') 
                       #print ('Plot arrows (dir)..') #(1, 141, 380, 1307) 
                       #for arr_idx_x in range(0,len(lons),10):
                       #    for arr_idx_y in range (0,len(lats),10):
                       #        if vals_du[arr_idx_y,arr_idx_x]!=0 and vals_dv[arr_idx_y,arr_idx_x]!=0:
                       #           arw = plt.arrow(lons[arr_idx_x],lats[arr_idx_y],vals_du[arr_idx_y,arr_idx_x]*3,vals_dv[arr_idx_y,arr_idx_x]*3,head_width=0.1 )

                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()


    ##################################
    # 2D VARIABLES :
    ##################################

    print ('I am working on 2D vars..')


    # ===============================
    # Loop on 2D VARS :
    # ===============================

    for idx_2d in range(0,len(field_2d_name)):
        var_2d=field_2d_name[idx_2d]
        var_2d_udm=field_2d_units[idx_2d]

        # ===============================
        # Loop on PERIOD/SEASONS :
        # ===============================
        for idx_dt in range(0,len(period_of_analysis)):
            dt_lab=period_of_analysis[idx_dt]

            print ('# 2D VAR [udm]: ',var_2d,' [',var_2d_udm,']')
            print ('# PERIOD: ',dt_lab)

            # Build the path/name of the nc file and open it 
            nc2open1=model_path+model_infileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname[0]+'.nc'
            nc2open2=model_path+model_infileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname[1]+'.nc'
            print ('Input files = ',nc2open1,nc2open2)
            model1 = NC.Dataset(nc2open1,'r')
            model2 = NC.Dataset(nc2open2,'r')
            #model1 = model1_file.mod_depth_var[vlev][:]
            #model2 = model_2file.mod_depth_var[vlev][:]

            # Read lat, lon and fild values 
            if grid == 'uv2t':
               lons = model1.variables['lon'][:]
               lats = model1.variables['lat'][:]
            else:
               lons = model1.variables['nav_lon'][:]
               lats = model1.variables['nav_lat'][:]
            vals = model1.variables[var_2d][:]-(model2.variables[var_2d][:]-0.25)
            vals_land=model1.variables[var_2d][:]
            lats_red=[]
            lons_red=[]
            vals_red=[]
            # mask land and fix max and/or min values 
            #thresh_inf = field_2d_inf[idx_2d]
            #thresh_sup = field_2d_sup[idx_2d]
            thresh = 0.0000
            mask = np.abs(vals_land) == thresh
            vals_ma = np.ma.masked_where(mask, vals)
            vals_ma_nan = np.ma.filled(vals_ma, -999.999)

            # Built the palette 
            vals_max=np.amax(vals_ma)
            vals_min=np.amin(vals_ma)
            #vals_max=10.0
            #vals_min=-10.0
            #vals_thresc=min(abs(vals_max),abs(vals_min))
            vals_thresc=(abs(vals_max)+abs(vals_min))/2
            #vals_thresc=0.25

            # Plot the map and save in the path/name
            plotname=model_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
            print ('Plot path/name: ',plotname)

            plt.figure(figsize=(15,8))
            plt.rc('font', size=12)
            # Plot Title
            if dt_lab == 'yr':
               plt.title ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_2d+' ['+var_2d_udm+']'+' - DEPTH: 0 [m] - DT: '+dates_label)
            else:
               plt.title ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_2d+' ['+var_2d_udm+']'+' - DEPTH: 0 [m] - DT: '+dt_lab+' ('+dates_label+')')
            # Read the coordinates for the plot 
            lon_0 = lons.mean()
            llcrnrlon = lons.min()
            urcrnrlon = lons.max()
            lat_0 = lats.mean()
            llcrnrlat = lats.min()
            urcrnrlat = lats.max()
            # Native vars case
            if grid != 'uv2t':
               # Create the map
               m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
               xi, yi = m(lons, lats)
               # Plot the frame to the map
               plt.rcParams["axes.linewidth"]  = 1.25
               # Contour of land with black line
               #land_value = [-999.999]
               # Plot the map
               #cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
               cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=vals_min,vmax=vals_max)
               # For sossheig + 25 cm:
               #vals_thresc=0.10 
               #cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
               text_min_x,text_min_y= m(-17,29.0)
               text_max_x,text_max_y= m(32,29.0)
               plt.text(text_max_x,text_max_y,'max='+str(vals_max), fontsize=12)
               plt.text(text_min_x,text_min_y,'min='+str(vals_min), fontsize=12)
               # Add the grid
               m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
               m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
               # Plot the legend and its label
               cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
               bar_label_string=var_2d+' ['+var_2d_udm+']'
               cbar.set_label(bar_label_string)
               # Land contour or color
               #land_value = [-999.999,vals_min-100]
               #m.contourf(xi,yi,np.squeeze(vals_ma_nan),land_value,colors='black')
               m.drawmapboundary(fill_color='gray')
               # Save and close 
               plt.savefig(plotname)
               plt.clf()
            # Derived var (i) case
            else:
                   # Plot the map
                   plt.xlim(llcrnrlon,urcrnrlon)
                   plt.ylim(llcrnrlat,urcrnrlat)
                   plt.xlabel('LON')
                   plt.ylabel('LAT')
                   print ('Plot colors (i)..')
                   cs = plt.pcolor(lons,lats,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
                   # Add the grid
                   plt.grid()
                   # Plot the legend and its label
                   land_value = -999.999
                   plt.colorbar(mappable=None,label='Currents'+' ['+var_3d_udm+']',orientation='horizontal',extend='max',aspect=60)
                   # Add land contour
                   if grid != 'uv2t':
                      contour = plt.contour(lons,lats,np.squeeze(vals_ma_nan),[land_value-500,land_value+500],colors='black')
                   else:
                      contour = plt.contour(lons,lats,np.squeeze(vals_ma_nan),[land_value-500,0.00],colors='black')
                   vals_max=np.amax(vals_ma)
                   vals_min=np.amin(vals_ma)
                   #text_min_x,text_min_y= m(-17,29.0)
                   #text_max_x,text_max_y= m(32,29.0)
                   text_max_x=32.0
                   text_max_y=29.0
                   text_min_x=-17.0
                   text_min_y=29.0
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,2)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)

                   # Save and close 
                   plt.savefig(plotname)
                   plt.clf()


