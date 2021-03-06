#
# Script for MAPS plot  
# You are supposed to provide a single nc file with the mean of the field on a single period and/or single file with field averaged on seasons. These can be built by map_extr.sh shell script.  
#
# by AC Goglio November 2019
#
#source activate mappyenv # activate into my virtual environment
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
from mpl_toolkits.basemap import Basemap # For plotting maps
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run paraeters:
# General run parameters:
#---------------------
grid = 'T' # Choose T, U or V grid
num_of_models = 2 # 1 for dataset maps or 2 for diffs between 2 datasets
dataset_name = 'diff' # '8' for Tides_8 run, 'ctrl' for the Control run or 'diff' for computing the diffs
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
#workdir_path = '/work/ag15419/tmp/map_ana_'+grid+dataset_name+'/'
workdir_path = '/work/ag15419/tmp/map_ana_'+grid+'diff'+'/'
#
# dates
#---------------------
# Choose start and end dates of the period 
inidate = '01/01/2016'
enddate = '31/12/2018'
# Choose periods for the analysis (yr means on the whole period), single file with seasonal mean should exists!
period_of_analysis=('yr','DJF','MAM','JJA','SON')
#--------------------
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
print ('Periods of analysis: ',period_of_analysis)
#
# field(s) to be plotted
#-----------------------
if grid == 'T':
   # T grid
   # 3D fields:
   field_3d_name=[] # #'votemper','vosaline'
   field_3d_units=['degC','PSU' ] # 'degC','PSU'
   field_3d_lev=[1.05,8,30,100,150,300,600,1000,2000] #,8,30,100,150,300,600,1000,2000] # 1.05,8,30,100,150,300,600,1000,2000
   # Fix thereshold for max and min map values for 3d fields or not
   T_threshold_3d=0
   field_3d_inf=[10,36]
   field_3d_sup=[25,42]
   # 2D fields:
   field_2d_name=['sossheig','sowaflup','soevapor','soprecip','sorunoff','soshfldo','sohefldo','solofldo','sosefldo','solafldo','somxl010']
   field_2d_units=['m','kg/m2/s','kg/m2/s','kg/m2/s','kg/m2/s','W/m2','W/m2','W/m2','W/m2','W/m2','m']
   field_2d_lev=0
   # Fix thereshold for max and min map values for 2d fields or not
   T_threshold_2d=0
   field_2d_inf=[10,36]
   field_2d_sup=[25,42]

elif grid == 'U':
     # U grid
     field_3d_name=['vozocrtx']
     field_3d_inf=[-0.4]
     field_3d_sup=[0.4]
     field_3d_units=['m/s']
     field_3d_lev=[1.05,8,30,100,150,300,600,1000,2000]
     #
     field_2d_name=['sozotaux']
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=['N/m2']
elif grid == 'V':
     # V grid
     field_3d_name=['vomecrty']
     field_3d_inf=[-0.3]
     field_3d_sup=[0.3]
     field_3d_units=['m/s']
     field_3d_lev=[1.05,8,30,100,150,300,600,1000,2000]
     #
     field_2d_name=['sometauy']
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=['N/m2']
elif grid == 'uv2t':
     # uv2t grid
     field_3d_name=['i']
     field_3d_inf=[0]
     field_3d_sup=[1]
     field_3d_units=['m/s']
     field_3d_lev=[1.05]
     #
     field_2d_name=[]
     field_2d_inf=[]
     field_2d_sup=[]
     field_2d_lev=0
     field_2d_units=[]
#
# SUB AREAS for currents analysis
# Gibraltar, Adriatic and Eastern Med Areas
Area_G_lon=[-7,9] # GB LONs 
Area_G_lat=[34,45] # GB LATs 
Area_A_lon=[11,23] # AD LONs 
Area_A_lat=[37,46] # AD LATs 
Area_EM_lon=[9,37] # EM LONs
Area_EM_lat=[30,41] # EM LATs
Area_ME_lon=[14,17] # ME LONs
Area_ME_lat=[36.5,39] # ME LATs
#
Area_name_ar=['Gibraltar','Adriatic','Eastern_Med','Messina']
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
      model_postname='mod_Control_run' #'mod_Tides8' or 'mod_Control_run'
      model_label='Control run' #'Tides_8 run' or 'Control run'
   # map3D_yr_%LEV%_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   print ('Model file templates: ',model_path,model_fileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_',model_postname,'.nc')
elif num_of_models == 2:
   print ('Model datasets DIFF MAPS..')
   model_path=workdir_path
   model_infileprename='map'
   model_fileprename='diff_map'
   model_postname=['mod_Tides8','mod_Control_run']
   model_label=['Tides_8 run','Control run']
   # map3D_yr_%LEV%_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   print ('Model file templates: ',model_path,model_infileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_','%model_postname%','.nc')

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
    
            # ===============================
            # Loop on PERIOD/SEASONS :
            # ===============================
            for idx_dt in range(0,len(period_of_analysis)):
                dt_lab=period_of_analysis[idx_dt]
    
                print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
                print ('# VERTICAL LEVEL: ',vlev)
                print ('# PERIOD: ',dt_lab)
    
                # Build the path/name of the nc file and open it 
                nc2open=model_path+model_fileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.nc'
                print ('Input file = ',nc2open)
                model = NC.Dataset(nc2open,'r')
    
                # Read lat, lon and fild values 
                if grid == 'uv2t':
                   lons = model.variables['lon'][:]
                   lats = model.variables['lat'][:]
                else: 
                   lons = model.variables['nav_lon'][:]
                   lats = model.variables['nav_lat'][:]
                vals = model.variables[var_3d][:]
                lats_red=[]
                lons_red=[]
                vals_red=[]
                # mask land and fix max and/or min values 
                thresh_inf = field_3d_inf[idx_3d]
                thresh_sup = field_3d_sup[idx_3d]
                thresh = 0.0000
                mask = np.abs(vals) == thresh 
                vals_ma = np.ma.masked_where(mask, vals)
                vals_ma_nan = np.ma.filled(vals_ma, -999.999)
                mask_rev= np.abs(vals) != thresh
                vals_ma_rev = np.ma.masked_where(mask_rev, vals)                
       
                # Plot the map and save in the path/name
                plotname=model_path+model_fileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.jpg'  
                print ('Plot path/name: ',plotname)
    
                plt.figure(figsize=(20,10))
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
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)

                # Derived var (i) case
                else:
                   # Plot the map
                   plt.xlim(llcrnrlon,urcrnrlon)
                   plt.ylim(llcrnrlat,urcrnrlat)
                   plt.xlabel('LON')
                   plt.ylabel('LAT')
                   cs = plt.pcolor(lons,lats,np.squeeze(vals_ma),cmap='jet',vmin=0.0,vmax=thresh_sup)
                   # Add the grid
                   plt.grid()
                   # Plot the legend and its label
                   land_value = 999.999
                   plt.colorbar(mappable=None,label='Currents'+' ['+var_3d_udm+']',orientation='horizontal',extend='max',aspect=60)
                   # Add land contour
                   contour = plt.contour(lons,lats,np.squeeze(vals_ma_nan),[land_value-500,land_value+100],colors='black')
                   vals_max=np.amax(vals_ma)
                   vals_min=np.amin(vals_ma)
                   text_min_x,text_min_y= m(-17,29.0)
                   text_max_x,text_max_y= m(32,29.0)
                   plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
                   plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)

                # Save and close 
                plt.savefig(plotname)
                plt.clf()
    
                # --------------------------------------------------
                # Plot the 3 map in different Med areas for currents
                if grid == 'U' or grid == 'V':

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
                       cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=Area_vmin_ar[idx_area],vmax=Area_vmax_ar[idx_area])
                       # Add the grid
                       m.drawparallels(np.arange(llcrnrlat,urcrnrlat, 1.), labels=[1,0,0,0], fontsize=10)
                       m.drawmeridians(np.arange(llcrnrlon,urcrnrlon, 2.), labels=[0,0,0,1], fontsize=10)
                       # Plot the legend and its label
                       cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
                       bar_label_string=var_3d+' ['+var_3d_udm+']'
                       cbar.set_label(bar_label_string)
                       vals_max=np.amax(vals_ma)
                       vals_min=np.amin(vals_ma)
                       #text_min_x,text_min_y= m(-17,29.0)
                       #text_max_x,text_max_y= m(32,29.0)
                       #plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
                       #plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)
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
            nc2open=model_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname+'.nc'
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
    
            plotname=model_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname+'.jpg'
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
               cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='jet',vmin=thresh_inf,extend='min')
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
            plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
            plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)
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

            # ===============================
            # Loop on PERIOD/SEASONS :
            # ===============================
            for idx_dt in range(0,len(period_of_analysis)):
                dt_lab=period_of_analysis[idx_dt]

                print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
                print ('# VERTICAL LEVEL: ',vlev)
                print ('# PERIOD: ',dt_lab)

                # Build the path/name of the nc file and open it 
                nc2open1=model_path+model_infileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'.nc'
                nc2open2=model_path+model_infileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[1]+'.nc'
                print ('Input files = ',nc2open1,nc2open2)
                model1 = NC.Dataset(nc2open1,'r')
                model2 = NC.Dataset(nc2open2,'r')

                # Read lat, lon and fild values 
                if grid == 'uv2t':
                   lons = model.variables['lon'][:]
                   lats = model.variables['lat'][:]
                else:
                   lons = model1.variables['nav_lon'][:]
                   lats = model2.variables['nav_lat'][:]
                vals = model1.variables[var_3d][:]-model2.variables[var_3d][:]
                lats_red=[]
                lons_red=[]
                vals_red=[]
                # mask land and fix max and/or min values 
                thresh_inf = field_3d_inf[idx_3d]
                thresh_sup = field_3d_sup[idx_3d]
                thresh = 0.0000
                #mask_sup = np.abs(vals) <= thresh_inf 
                #mask_inf = np.abs(vals) >= thresh_sup
                #vals_ma = np.ma.masked_where(mask_sup, vals)
                #vals_ma = np.ma.masked_where(mask_inf, vals_ma)
                mask = np.abs(vals) == thresh
                vals_ma = np.ma.masked_where(mask, vals)
                vals_ma_nan = np.ma.filled(vals_ma, -999.999)

                # Built the palette 
                vals_max=np.amax(vals_ma)
                vals_min=np.amin(vals_ma)
                if var_3d == 'vosaline':
                   vals_thresc=min(abs(vals_max),abs(vals_min))
                elif str(vlev) == '8' or str(vlev) == '30'or str(vlev) == '100':
                     vals_thresc=min(abs(vals_max),abs(vals_min))
                else:
                   vals_thresc=(abs(vals_max)+abs(vals_min))/2

                # Plot the map and save in the path/name

                plotname=model_path+model_fileprename+'3D_'+dt_lab+'_'+str(vlev)+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
                print ('Plot path/name: ',plotname)

                plt.figure(figsize=(20,10))
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
                # Create the map
                m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
                xi, yi = m(lons, lats)
                # Plot the frame to the map
                plt.rcParams["axes.linewidth"]  = 1.25
                # Contour of land with black line
                land_value = [0.0000]
                #contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
                # Plot the map
                #cmap = plt.cm.bwr
                #cmap.set_bad('black',1.)
                cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
                #cs = m.contourf(xi,yi,np.squeeze(vals))
                text_min_x,text_min_y= m(-17,29.0)
                text_max_x,text_max_y= m(32,29.0)
                plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
                plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)
                # Add the grid
                m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
                m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
                #m.drawcoastlines()
                #m.fillcontinents(color='white')
                # Plot the legend and its label
                cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
                bar_label_string=var_3d+' ['+var_3d_udm+']'
                cbar.set_label(bar_label_string)
                # Land contour
                land_value = [-999.999,vals_min-100]
                m.contourf(xi,yi,np.squeeze(vals_ma_nan),land_value,colors='black')
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

            # Read lat, lon and fild values 
            if grid == 'uv2t':
               lons = model.variables['lon'][:]
               lats = model.variables['lat'][:]
            else:
               lons = model1.variables['nav_lon'][:]
               lats = model1.variables['nav_lat'][:]
            vals = model1.variables[var_2d][:]-model2.variables[var_2d][:]
            lats_red=[]
            lons_red=[]
            vals_red=[]
            # mask land and fix max and/or min values 
            thresh_inf = field_2d_inf[idx_2d]
            thresh_sup = field_2d_sup[idx_2d]
            thresh = 0.0000
            #mask_sup = np.abs(vals) <= thresh_inf 
            #mask_inf = np.abs(vals) >= thresh_sup
            #vals_ma = np.ma.masked_where(mask_sup, vals)
            #vals_ma = np.ma.masked_where(mask_inf, vals_ma)
            mask = np.abs(vals) == thresh
            vals_ma = np.ma.masked_where(mask, vals)
            vals_ma_nan = np.ma.filled(vals_ma, -999.999)

            # Built the palette 
            vals_max=np.amax(vals_ma)
            vals_min=np.amin(vals_ma)
            
            vals_thresc=min(abs(vals_max),abs(vals_min))
            #vals_thresc=(abs(vals_max)+abs(vals_min))/2


            # Plot the map and save in the path/name
            plotname=model_path+model_fileprename+'2D_'+dt_lab+'_0_'+var_2d+'_'+dates_label+'_'+model_postname+'.jpg'
            print ('Plot path/name: ',plotname)

            plt.figure(figsize=(20,10))
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
            # Create the map
            m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
            xi, yi = m(lons, lats)
            # Plot the frame to the map
            plt.rcParams["axes.linewidth"]  = 1.25
            # Contour of land with black line
            land_value = [0.0000]
            #contour = plt.contour(xi,yi,np.squeeze(vals),land_value, colors='black')
            # Plot the map
            cs = m.pcolor(xi,yi,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
            #cs = m.contourf(xi,yi,np.squeeze(vals))
            text_min_x,text_min_y= m(-17,29.0)
            text_max_x,text_max_y= m(32,29.0)
            plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=12)
            plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,1)), fontsize=12)
            # Add the grid
            m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
            m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
            #m.drawcoastlines()
            #m.fillcontinents(color='white')
            # Plot the legend and its label
            cbar = m.colorbar(cs, location='bottom', pad="10%",extend='both')
            bar_label_string=var_2d+' ['+var_2d_udm+']'
            cbar.set_label(bar_label_string)
            # Land contour
            land_value = [-999.999,vals_min-1]
            m.contourf(xi,yi,np.squeeze(vals_ma_nan),land_value,colors='black')
            # Save and close 
            plt.savefig(plotname)
            plt.clf()

