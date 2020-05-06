#
# Script for SECTION plot  
# You are supposed to provide a single nc file with the mean of the field on a single period and/or single file with field averaged on seasons. These can be built by map_extr.sh shell script.  
#
# by AC Goglio November 2019
#
#source activate mappyenv # activate my virtual environment
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
from matplotlib.gridspec import GridSpec # For plot division in sublots
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
grid = 'T' # If grid == W remember to modify the digits number and the position of the plt.text, reduce the coo range  
num_of_models = 1 # 1 for dataset maps or 2 for diffs between 2 datasets
dataset_name = 'ctrl' # '8' for Tides_8 run, 'ctrl' for the Control run or 'diff' for computing the diffs
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
workdir_path = '/work/ag15419/tmp/sect_ana_'+grid+'_'+dataset_name
#workdir_path = '/work/ag15419/tmp/sect_ana_'+grid+dataset_name+'_me/'
#workdir_path = '/work/ag15419/tmp/sect_ana_Tdiff_me'+'/'
#workdir_path = '/work/ag15419/tmp/sect_ana_'+grid+'diff'+'/'
#
# dates
#---------------------
# Choose start and end dates of the period 
inidate = '01/01/2016'
enddate = '31/12/2016'
# Choose periods for the analysis (yr means on the whole period), single file with seasonal mean should exists!
period_of_analysis=('yr','DJF','MAM','JJA','SON')
#--------------------
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
print ('Periods of analysis: ',period_of_analysis)
#
# Section infos
#--------------------
lat_range='ALL' # 'ALL' or 2elements array
lon_range='ALL'  # 'ALL' or 2elements array, e.g. [-8.000,-4.000] 
sect_label='GB' # GB,SI,ME,GBv
sect_name='Messina' #'Gibraltar' 'Sicily','Messina'
depth_first_layer=-300
depth_inf=-1500
mod_lev_inf=-85 # -85 for GB, ME, SI or -65 for GBv
print ('Section name',sect_name,' (',sect_label,')')
#-----------------------
# field(s) to be plotted
if grid == 'T':
   # T grid
   field_3d_name=['votemper','vosaline'] # 'votemper','vosaline'
   field_3d_inf=[12,35]
   field_3d_sup=[20,39]
   field_3d_units=['degC','PSU' ] # 'degC','PSU'
   mod_depth_var='deptht'
   #
elif grid == 'U':
     # U grid
     field_3d_name=['vozocrtx']
     field_3d_inf=[-0.6]
     field_3d_sup=[0.6]
     field_3d_units=['m/s']
     mod_depth_var='depthu'
     #
elif grid == 'V':
     # V grid
     field_3d_name=['vomecrty']
     field_3d_inf=[-0.4]
     field_3d_sup=[0.4]
     field_3d_units=['m/s']
     mod_depth_var='depthv'
     #
elif grid == 'W':
     # W grid
     field_3d_name=['vovecrtz']
     field_3d_inf=[-0.01]
     field_3d_sup=[0.01]
     field_3d_units=['m/s']
     mod_depth_var='depthw'
elif grid == 'uv2t':
     # W grid
     field_3d_name=['i']
     field_3d_inf=[0.0]
     field_3d_sup=[1.0]
     field_3d_units=['m/s']
     mod_depth_var='depth'
#
# For the text in the plot
if sect_label == 'GB':
   if lon_range[1] == -4.0:
        text_pos=[-8,-4.5]
   else:
        text_pos=[-10,-1]
elif sect_label == 'SI':
     text_pos=[33,41.5]
elif sect_label == 'ME':
     text_pos=[38,38.45]
elif sect_label == 'GBv':
     text_pos=[35.90,36.00]
     #
if grid == 'W':
   text_dig=3
elif sect_label == 'ME':
   text_dig=2
elif sect_label == 'SI':
   text_dig=2
else:
   text_dig=1

#
# INPUTS
# Path and name of inputs datasets
# Currently the extraction of nc is done externally by another script
#
# MODEL DATASETS
# Currently this scripts plots maps for just 1 dataset (can be improoved..)
# Model 1st db file template (at the moment just one model db is possible): %model1obs_path%/%model1obs_prename%_stn_%model1obs_postname%.nc
if num_of_models == 1:
   print ('Model dataset..')
   if dataset_name == '8':
      model_path=workdir_path
      model_fileprename='sect'
      model_postname='mod_Tides8' #'mod_Tides8' or 'mod_Control_run'
      model_label='Tides_8 run' #'Tides_8 run' or 'Control run'
   elif dataset_name == 'ctrl':
      model_path=workdir_path
      model_fileprename='sect'
      model_postname='mod_Control_run' #'mod_Tides8' or 'mod_Control_run'
      model_label='Control run' #'Tides_8 run' or 'Control run'
   # sect_%SECT_LABEL%_yr_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   print ('Model file templates: ',model_path,model_fileprename,'%section_label%','%periods_of_analysis%','%field%',dates_label,'_',model_postname,'.nc')
elif num_of_models == 2:
   print ('Model datasets DIFF SECTION..')
   model_path=workdir_path
   model_infileprename='sect'
   model_fileprename='diff_sect'
   model_postname=['mod_Tides8','mod_Control_run']
   model_label=['Tides_8 run','Control run']
   thresh_f=-0.0 # To have a narrow palette (for U = -0.05; for T = 0.0 )
   # sect_%SECT_LABEL%_yr_%FIELD%_%ANA_STARTDATE%_%ANA_ENDDATE%_mod_%INDATASET%.nc
   print ('Model file templates: ',model_path,model_fileprename,'%section_label%','%periods_of_analysis%','%field%',dates_label,'_','%model_postname%','.nc')

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
print ('Fields: ',field_3d_name)
print ('Periods (global+season): ',dates_label,period_of_analysis)

if num_of_models == 1:

    print ('I am working on the following model dataset: ',model_label)
    print ('I am extraction the following section: ',sect_name)
    # ===============================
    # Loop on VARS :
    # ===============================
    
    for idx_3d in range(0,len(field_3d_name)):
        var_3d=field_3d_name[idx_3d]
        var_3d_udm=field_3d_units[idx_3d]
    
        # ===============================
        # Loop on PERIOD/SEASONS :
        # ===============================
        for idx_dt in range(0,len(period_of_analysis)):
            dt_lab=period_of_analysis[idx_dt]
  
            print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
            print ('# PERIOD: ',dt_lab)

            # Build the path/name of the nc file and open it 
            nc2open=model_path+model_fileprename+'_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.nc'
            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon and fild values 
            lons = model.variables['nav_lon'][:]
            lats = model.variables['nav_lat'][:]
            depths = model.variables[mod_depth_var][:]
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

            # Plot the map and save in the path/name
            if lat_range == 'ALL' and lon_range == 'ALL':
               plotname=model_path+model_fileprename+'_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.jpg'  
            else:
               plotname=model_path+model_fileprename+'_'+'sub_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname+'.jpg'
            print ('Plot path/name: ',plotname)

            fig = plt.figure(figsize=(12,8))
            plt.rc('font', size=14)
            # Plot Title
            if dt_lab == 'yr':
               fig.suptitle ('Model '+model_label+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
               print ('Model '+model_label+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
            else:
               fig.suptitle ('Model '+model_label+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')
               print ('Model '+model_label+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')            

            # Read the coordinates for the plot 
            # Longitude is constant 
            if lons.min() == lons.max():       
               print ('Longitude is constant')
               coo_0 = lats.mean()
               if lat_range == 'ALL':
                  llcrnrcoo = lats.min()
                  urcrnrcoo = lats.max()
               else:
                  llcrnrcoo = lat_range[0]
                  urcrnrcoo = lat_range[1]
               coos = lats
               print(lons.min(),llcrnrcoo,urcrnrcoo)
               xlabel_string='LAT '+'('+'LON = '+str(lons.min())+')'
               #
            # Latitude is constant
            elif lats.min() == lats.max():
                 print ('Latitude is constant')
                 coo_0 = lons.mean()
                 if lon_range == 'ALL':
                    llcrnrcoo = lons.min()
                    urcrnrcoo = lons.max()
                 else:
                    llcrnrcoo = lon_range[0]
                    urcrnrcoo = lon_range[1]
                 coos = lons
                 xlabel_string='LON '+'('+'LAT = '+str(lats.min())+')'
                 #
            depth_0 = depths.mean()
            llcrnrdepth = depths.min()
            urcrnrdepth = depths.max()
                 
            # Section of first layers
            ax=plt.subplot(1,1,1)

            # Plot the frame to the map
            plt.rcParams["axes.linewidth"]  = 1.25
            # Plot the map
            print ('Plotting..')
            ###plt.ylabel('Depth [m]')
            plt.ylabel('Mod lev [#]')
            plt.xlabel(xlabel_string)
            ###plt.ylim(depth_first_layer,0.0)
            plt.ylim(mod_lev_inf,0.0)
            depths_inv=-1.0*depths
            mod_lev=np.arange(0,len(depths),1)
            plt.xlim(llcrnrcoo,urcrnrcoo)
            # comp with Cucco in Messina strait 
            if var_3d == 'votemper':
               vmin=14.8
               vmax=21.2
            elif var_3d == 'vosaline':
                 vmin=38.0
                 vmax=39.0
            cs = plt.pcolor(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma),cmap='gist_rainbow',vmin=vmin,vmax=vmax)

            # Built the palette  (NON Cucco)
            vals_max=np.amax(vals_ma)
            vals_min=np.amin(vals_ma)
            #cs = plt.pcolor(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma),cmap='jet')
            #
            # Plot the legend and its label
            plt.colorbar(mappable=None,label=var_3d+' ['+var_3d_udm+']',orientation='horizontal',extend='both')
            # Contour of land with black line
            land_value=-999.999
            ###contour = plt.contourf(np.squeeze(coos),depths_inv,np.squeeze(vals_ma_nan),[land_value,land_value+1],colors='gray')
            contour = plt.contourf(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma_nan),[land_value-500,thresh_inf-1],colors='gray')
            # Add second y axes with depth
            ax.set_yticklabels(['ff','80','70','60','50','40','30','20','10','0'])
            plt.grid()
            ax2 = ax.twinx()
            ax2.set_ylim(mod_lev_inf,0.0) 
            ax2.set_ylabel('Depth [m]')
            ax2.set_yticklabels(['ff','-1229 m','-874 m','-602 m','-399 m','-250 m','-145 m','-73 m','-26 m'])
            # Add max and min palette values
            plt.text(text_pos[1],-100,'max='+str(round(vals_max,text_dig)), fontsize=14) 
            plt.text(text_pos[0],-100,'min='+str(round(vals_min,text_dig)), fontsize=14)
            # Add the grid
            #plt.grid()

#            # Section of all the layers
#            plt.subplot(2,1,2)
#
#            # Plot the frame to the map
#            plt.rcParams["axes.linewidth"]  = 1.25
#            # Contour of land with black line
#            land_value=0.0
#            # Plot the map
#            print ('Plotting..')
#            plt.xlabel(xlabel_string)
#            plt.ylabel('Depth [m]')
#            plt.ylim(depth_inf,0.0)
#            depths_inv=-1.0*depths
#            plt.xlim(llcrnrcoo,urcrnrcoo)
#            # Built the palette 
#            vals_max=np.amax(vals_ma)
#            vals_min=np.amin(vals_ma)
#            #palette = plt.cm.jet
#            #palette.set_bad (color='black')
#            cs = plt.pcolor(coos,depths_inv,np.squeeze(vals_ma),cmap='jet',vmin=thresh_inf,vmax=thresh_sup) 
#            # Plot the legend and its label
#            plt.colorbar(label=var_3d+' ['+var_3d_udm+']')
#            # Contour of land with black line
#            land_value=999.999
#            contour = plt.contourf(np.squeeze(coos),depths_inv,np.squeeze(vals_ma_nan),[land_value,land_value+1],colors='gray')
#            # Add max and min palette values
#            plt.text(text_pos,0.0,'max='+str(round(vals_max,text_dig)), fontsize=12)
#            plt.text(text_pos,depth_inf,'min='+str(round(vals_min,text_dig)), fontsize=12)
#            # Add the grid
#            plt.grid()

            # Save and close 
            plt.savefig(plotname)
            plt.clf()






##########################################
########## DIFF BETWEEN MODELS ###########

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
   print ('Diff case!')

   # ===============================
   # Loop on 3D VARS :
   # ===============================
   
   for idx_3d in range(0,len(field_3d_name)):
       var_3d=field_3d_name[idx_3d]
       var_3d_udm=field_3d_units[idx_3d]
   
       # ===============================
       # Loop on PERIOD/SEASONS :
       # ===============================
       for idx_dt in range(0,len(period_of_analysis)):
               dt_lab=period_of_analysis[idx_dt]
   
               print ('# 3D VAR [udm]: ',var_3d,' [',var_3d_udm,']')
               print ('# PERIOD: ',dt_lab)
   
               # Build the path/name of the nc file and open it 
               nc2open1=model_path+model_infileprename+'_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'.nc'
               nc2open2=model_path+model_infileprename+'_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname[1]+'.nc'
               print ('Input files = ',nc2open1,nc2open2)
               model1 = NC.Dataset(nc2open1,'r')
               model2 = NC.Dataset(nc2open2,'r')
   
               # Read lat, lon and fild values 
               if grid == 'uv2t':
                  lons = model1.variables['lon'][:]
                  lats = model1.variables['lat'][:]
               else:
                  lons = model1.variables['nav_lon'][:]
                  lats = model1.variables['nav_lat'][:]

               depths = model1.variables[mod_depth_var][:]
               vals = model1.variables[var_3d][:]-model2.variables[var_3d][:]
               vals_land = model1.variables[var_3d][:]

               # mask land and fix max and/or min values 
               thresh_inf = field_3d_inf[idx_3d]
               thresh_sup = field_3d_sup[idx_3d]
               thresh = 0.0000
               mask = np.abs(vals) == thresh
               vals_ma = np.ma.masked_where(mask, vals_land)
               vals_ma_nan = np.ma.filled(vals_ma, -999.999)
               valsdiff_ma = np.ma.masked_where(mask, vals ) 

               # Plot the map and save in the path/name
               if lat_range == 'ALL' and lon_range == 'ALL':
                  plotname=model_path+model_fileprename+'_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
               else:
                  plotname=model_path+model_fileprename+'_sub_'+sect_label+'_'+dt_lab+'_'+var_3d+'_'+dates_label+'_'+model_postname[0]+'_'+model_postname[1]+'.jpg'
               print ('Plot path/name: ',plotname)

               # Built the palette 
               vals_max=np.amax(vals_ma)
               vals_min=np.amin(vals_ma)
               if var_3d == 'votemper':
                  #vals_thresc=(abs(vals_max)+abs(vals_min))/2
                  vals_thresc=max(abs(vals_max),abs(vals_min))
               elif var_3d == 'vovecrtz':
                  #vals_thresc=(abs(vals_max)+abs(vals_min))/2
                  vals_thresc=0.006
               else:
                  #vals_thresc=min(abs(vals_max),abs(vals_min))+thresh_f
                  vals_thresc=max(abs(vals_max),abs(vals_min))
                  #vals_thresc=(abs(vals_max)+abs(vals_min))/2
               # TO BE REMOVED:
               if var_3d == 'votemper':
                  vals_thresc=0.5
               elif var_3d == 'vosaline':
                  vals_thresc=0.3

               # Read the coordinates for the plot 
               # Longitude is constant 
               if lons.min() == lons.max():
                  print ('Longitude is constant')
                  coo_0 = lats.mean()
                  if lat_range == 'ALL':
                     llcrnrcoo = lats.min()
                     urcrnrcoo = lats.max()
                  else:
                     llcrnrcoo = lat_range[0]
                     urcrnrcoo = lat_range[1]
                  coos = lats
                  xlabel_string='LAT '+'('+'LON = '+str(lons.min())+')'
                  #
               # Latitude is constant
               elif lats.min() == lats.max():
                    print ('Latitude is constant')
                    coo_0 = lons.mean()
                    if lon_range == 'ALL':
                       llcrnrcoo = lons.min()
                       urcrnrcoo = lons.max()
                    else:
                       llcrnrcoo = lon_range[0]
                       urcrnrcoo = lon_range[1]
                    coos = lons
                    xlabel_string='LON '+'('+'LAT = '+str(lats.min())+')'
                    #
               depth_0 = depths.mean()
               llcrnrdepth = depths.min()
               urcrnrdepth = depths.max()

               # Plot the map
               print ('Plotting..')
               # Colormap analysis (not arrow)
               if sect_label != 'GBv':
                  fig = plt.figure(figsize=(12,12))
                  plt.rc('font', size=14)
                  plt.rcParams["axes.linewidth"]  = 1.25
                  # Plot Title
                  if dt_lab == 'yr':
                     fig.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
                     print ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
                  else:
                     fig.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')
                     print ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')
                  # Section of first layers
                  ax=plt.subplot(1,1,1)
                  ###plt.ylabel('Depth [m]')
                  plt.ylabel('Mod lev [#]')
                  ###plt.ylim(depth_first_layer,0.0)
                  plt.ylim(mod_lev_inf,0.0)
                  depths_inv=-1.0*depths
                  mod_lev=np.arange(0,len(depths),1)
                  plt.rcParams["axes.linewidth"]  = 1.25
                  plt.xlabel(xlabel_string)
                  plt.xlim(llcrnrcoo,urcrnrcoo)
                  # Plot
                  if grid != 'uv2t':
                     # comp with Cucco in Messina strait 
                     if var_3d == 'votemper':
                        vmin=14.8
                        vmax=21.2
                     elif var_3d == 'vosaline':
                        vmin=38.0
                        vmax=39.0
                     cs = plt.pcolor(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma),cmap='gist_rainbow',vmin=vmin,vmax=vmax)
                     plt.colorbar(mappable=None,label=var_3d+' ['+var_3d_udm+']',orientation='horizontal',extend='both')
                     #
                     #cs = plt.pcolor(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
                     #plt.colorbar(mappable=None,label=var_3d+' ['+var_3d_udm+']',orientation='horizontal',extend='both')
                  elif grid == 'uv2t':
                     cs = plt.pcolor(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma),cmap='Reds',vmin=0,vmax=vals_thresc)
                     plt.colorbar(mappable=None,label=var_3d+' ['+var_3d_udm+']',orientation='horizontal',extend='max')
                  # Contour of land with black line
                  land_value=-999.999
                  ###contour = plt.contourf(np.squeeze(coos),depths_inv,np.squeeze(vals_ma_nan),[land_value-500,land_value+100],colors='gray')
                  #contour = plt.contourf(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma_nan),[land_value-500,-vals_thresc-10],colors='gray')
                  contourf = plt.contourf(np.squeeze(coos),-mod_lev,np.squeeze(vals_ma_nan),[land_value-300,land_value+300],colors='gray')
                  # Add second y axes with depth
                  ax.set_yticklabels(['ff','80','70','60','50','40','30','20','10','0'])
                  plt.grid()
                  ax2 = ax.twinx()
                  ax2.set_ylim(mod_lev_inf,0.0)
                  ax2.set_ylabel('Depth [m]')
                  ax2.set_yticklabels(['ff','-1229 m','-874 m','-602 m','-399 m','-250 m','-145 m','-73 m','-26 m'])
                  # Add max and min palette values
                  plt.text(text_pos[1],-100,'max='+str(round(vals_max,text_dig)), fontsize=14)
                  plt.text(text_pos[0],-100,'min='+str(round(vals_min,text_dig)), fontsize=14)

               else:
               # ARROW case in GBv
                  fig = plt.figure(figsize=(25,10))
                  plt.rc('font', size=14)
                 # Plot Title
                  if dt_lab == 'yr':
                     fig.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
                     print ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dates_label)
                  else:
                     fig.suptitle ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')
                     print ('DIFF '+model_label[0]+'-'+model_label[1]+' MEAN: '+var_3d+' ['+var_3d_udm+']'+' - SECTION: '+sect_name+' - DT: '+dt_lab+' ('+dates_label+')')
                  mod_lev_inf=-65
                  plt.ylabel('Mod lev [#]')
                  plt.ylim(mod_lev_inf,1)
                  depths_inv=-1.0*depths
                  mod_lev=np.arange(0,len(depths),1)
                  #
                  arrow_bases=np.squeeze(coos)
                  arrow_mean=np.zeros((141))
                  arrow_num=np.zeros((141))
                  arrow_net=np.zeros((141))
                  for arrow_base in range(0,len(arrow_bases)+1):
                      ax=plt.subplot(1,5,arrow_base+1)
                      plt.grid()
                      if arrow_base == 0:
                         plt.ylabel('Mod lev [#]')
                      plt.ylim(mod_lev_inf,1)
                      depths_inv=-1.0*depths
                      mod_lev=np.arange(0,len(depths),1)
                      plt.title('LON/LAT = ['+str(lons.min())+';'+str(arrow_base)+']')
                      plt.xlabel('Zonal currents diff [m/s]')
                      # Divide single points from mean and avg computation
                      if arrow_base != len(arrow_bases):
                         arrow_dim=np.squeeze(vals[:,:,arrow_base,:])
                         #mask_arrow=0.0
                         #arrow_dim_masked=np.ma.masked_where(mask_arrow,arrow_dim)
                         #print (arrow_dim_masked[-1])
                         arrow_mean=arrow_mean+arrow_dim
                         arrow_net=arrow_net+arrow_dim
                         #
                         plt.title('LON/LAT = ['+str(lons.min())+';'+str(arrow_bases[arrow_base])+']')
                         #plt.xlim(-np.amax(arrow_dim)-0.1,np.amax(arrow_dim)+0.1)
                         plt.xlim(-0.3,0.3)
                         for arrow_idx in range (0,len(arrow_dim)):
                             plt.arrow(0,-mod_lev[arrow_idx],arrow_dim[arrow_idx],0,head_width=1, head_length=0.01)
                             if arrow_dim[arrow_idx] != 0:
                                arrow_num[arrow_idx]=arrow_num[arrow_idx]+1
                      else:
                         arrow_mean=arrow_mean/(arrow_num) 
                         arrow_nan = np.ma.filled(arrow_mean,0.0)
                         #plt.title('NET')
                         plt.title('MEAN')
                         #plt.xlim(-np.amax(arrow_nan)-0.1,np.amax(arrow_nan)+0.1)
                         #plt.xlim(-0.5,0.5)
                         plt.xlim(-0.3,0.3)
                         for arrow_idx in range (0,len(arrow_nan)):
                             #plt.arrow(0,-mod_lev[arrow_idx],arrow_net[arrow_idx],0,head_width=1, head_length=0.01,fc='Red')
                             plt.arrow(0,-mod_lev[arrow_idx],arrow_nan[arrow_idx],0,head_width=1, head_length=0.01,fc='Red')
                         #
                         # Add second y axes with depth
                         if arrow_base == 0:
                            ax.set_yticklabels(['ff','60','50','40','30','20','10','0'])
                         ax2 = ax.twinx()
                         ax2.set_ylim(mod_lev_inf,1)
                         ax2.set_ylabel('Depth [m]')
                         ax2.set_yticklabels(['ff','-602 m','-399 m','-250 m','-145 m','-73 m','-26 m','0 m'])


#               # Section of all the layers
#               plt.subplot(2,1,2)
#               # Plot the frame to the map
#               plt.rcParams["axes.linewidth"]  = 1.25
#               # Plot the map
#               print ('Plotting..')
#               plt.xlabel(xlabel_string)
#               plt.ylabel('Depth [m]')
#               plt.ylim(depth_inf,0.0)
#               plt.xlim(llcrnrcoo,urcrnrcoo)
#               depths_inv=-1.0*depths
#               # Plot
#               cs = plt.pcolor(coos,depths_inv,np.squeeze(vals_ma),cmap='bwr',vmin=-vals_thresc,vmax=vals_thresc)
#               plt.colorbar(label=var_3d+' ['+var_3d_udm+']')
#               # Contour of land with black line
#               land_value=999.999
#               contour = plt.contourf(np.squeeze(coos),depths_inv,np.squeeze(vals_ma_nan),[land_value,land_value+1], colors='gray')
#               # Add max and min values infos
#               plt.text(text_pos,0.0,'max='+str(round(vals_max,text_dig)), fontsize=12)
#               plt.text(text_pos,depth_inf,'min='+str(round(vals_min,text_dig)), fontsize=12)
#               # Add the grid
#               plt.grid()

               # Save and close 
               plt.savefig(plotname)
               plt.clf()

