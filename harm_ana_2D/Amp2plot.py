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
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run paraeters:
# General run parameters:
#---------------------
# work dir path (WARNING: file in the directory will be removed..)
workdir_path = '/work/ag15419/tmp/HA_maps_andTG'+'/'
#
# dates
#---------------------
# Choose start and end dates of the period 
inidate = '01/07/2017'
enddate = '31/12/2017'
#--------------------
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
#
# field(s) to be plotted
#-----------------------
mod_depth_var='deptht'
# 2D fields:
field_2d_name=['M2_Amp'] #,'K1_Amp','O1_Amp','S2_Amp','P1_Amp','N2_Amp','Q1_Amp','K2_Amp']
P_field_2d_name=['M2_Pha','K1_Pha','O1_Pha','S2_Pha','P1_Pha','N2_Pha','Q1_Pha','K2_Pha']
field_2d_units=['cm','cm','cm','cm','cm','cm','cm','cm']
P_field_2d_units=['deg','deg','deg','deg','deg','deg','deg','deg']
field_2d_lev=0
# Fix thereshold for max map values for 2d fields or not
thsup_flag=1
#field_2d_sup=[50,25,10,25,5,5,2,3] # Valentina's values
field_2d_sup=[25,15,5,15,5,5,2,3]
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


# INPUTS
# Path and name of inputs datasets
# Currently the extraction of nc is done externally by another script
#
# MODEL DATASETS
# Currently this scripts plots maps for just 1 dataset (can be improoved..)
# Model 1st db file template (at the moment just one model db is possible): %model1obs_path%/%model1obs_prename%_stn_%model1obs_postname%.nc
model_path='/work/ag15419/arc_pp/tides/h_arc/map_ha/'
model_fileprename='amppha'
model_postname='mod_Tides8' #'mod_Tides8' or 'mod_Control_run'
model_label='Tides_8 run' #'Tides_8 run' or 'Control run'
print ('Model file templates: ',model_path,model_fileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_',model_postname,'.nc')

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
print ('Fields 2D: ',field_2d_name)

print ('I am working on the following model dataset: ',model_label)
    
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
    
            P_var_2d=P_field_2d_name[idx_2d]
            P_var_2d_udm=P_field_2d_units[idx_2d]

            print ('# 2D VAR [udm]: ',var_2d,' [',var_2d_udm,']')
    
            # Build the path/name of the nc file and open it 
            nc2open=model_path+model_fileprename+'2D_'+'0_'+'sossheig'+'_'+dates_label+'_'+model_postname+'.nc'
            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon and fild values 
            lons = model.variables['nav_lon'][:]
            lats = model.variables['nav_lat'][:]
            vals = model.variables[var_2d][:]*100 # Want cm not meters!
            P_vals = model.variables[P_var_2d][:]          

            lats_red=[]
            lons_red=[]
            vals_red=[]
            # mask land and fix max and/or min values 
            #thresh_inf = field_2d_inf[idx_2d]
            #thresh_sup = field_2d_sup[idx_2d]
            thresh = 0.0000
            #mask_sup = np.abs(vals) <= thresh_inf 
            #mask_inf = np.abs(vals) >= thresh_sup
            #vals_ma = np.ma.masked_where(mask_sup, vals)
            #vals_ma = np.ma.masked_where(mask_inf, vals_ma)
            mask = np.abs(vals) == thresh
            vals_ma = np.ma.masked_where(mask, vals)
            
            thresh2 = -199    
            mask2 = np.abs(vals) == thresh2
            vals_ma2 = np.ma.masked_where(mask2, vals)

            # Plot the map and save in the path/name
    
            plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Amp.jpg'
            print ('Plot path/name: ',plotname)
    
            plt.figure(figsize=(20,10))
            plt.rc('font', size=12)
            # Plot Title
            plt.title (var_2d+'litude ['+var_2d_udm+'] - Harmonic Analysis '+dates_label)
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
            # Plot the map and add the amphidromes
            if thsup_flag == 1:
               cs = m.pcolor(xi,yi,np.abs(np.squeeze(vals)),cmap='jet',vmax=field_2d_sup[idx_2d])
               # Valentina's phase analysis:
               amphid_p = plt.contour(xi,yi,np.squeeze(P_vals),[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360],colors='white')
               amphid_n = plt.contour(xi,yi,np.squeeze(P_vals),[-360,-340,-320,-300,-280,-260,-240,-220,-200,-180,-160,-140,-120,-100,-80,-60,-40,-20],colors='white',linestyles='dashed')
               #amphid = plt.contour(xi,yi,np.squeeze(P_vals),360, colors='white')
            else:
               cs = m.pcolor(xi,yi,np.abs(np.squeeze(vals)),cmap='jet')
               amphid = plt.contour(xi,yi,np.squeeze(P_vals),360, colors='white') 
            # Contour of land with black line 
            #land_value = [0.00000,10.0,20.0,30.0,40.0,50.0] # Valentina's contours
            land_value = 0.00000
            contour = plt.contour(xi,yi,np.abs(np.squeeze(vals_ma)),land_value, colors='black') #colors='black'
            # Add the grid
            m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
            m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
            # Plot the legend and its label
            cbar = m.colorbar(cs, location='bottom', pad="10%",extend='max')
            bar_label_string=var_2d+'litude'+' ['+var_2d_udm+']'
            cbar.set_label(bar_label_string)
            vals_max=np.amax(abs(vals_ma))
            vals_min=0
            #text_min_x,text_min_y= m(-17,29.0)
            text_max_x,text_max_y= m(32,29.0)
            plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1))+var_2d_udm, fontsize=12)
            #plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=12)
            #
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

#            # Atlantic Box
#            xp, yp = m(-1.51364,43.52696)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG AngletConvergentTG 
#            #xp, yp = m(-1.16350,44.66500)
#            #plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG ArcachonEyracTG 
#            xp, yp = m(-1.51540,43.52730)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG BayonneBoucauTG 
#            #xp, yp = m(-1.47315,43.48183)
#            #plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG BayonnePontBlancTG 
#            xp, yp = m(-1.47941,43.49758)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG BayonneQuaiDeLessepsTG 
#            xp, yp = m(-3.05,43.357)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG BilbaoTG 
#            xp, yp = m(-6.34,36.8)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#ff1493') # TG BonanzaTG 
#            xp, yp = m(-1.66271,43.38916)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG CiboureTG 
#            xp, yp = m(-8.389,43.357)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG CorunaTG 
#            xp, yp = m(-8.2485,43.4762)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG Ferrol2TG 
#            xp, yp = m(-8.326,43.463)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG FerrolTG 
#            xp, yp = m(-5.698,43.558)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG GijonTG 
#            xp, yp = m(-6.834,37.132)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#ff1493') # TG HuelvaTG 
#            xp, yp = m(-8.5301,43.3465)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG LangosteiraTG 
#            xp, yp = m(-8.7045,41.1867)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG LeixoesTG 
#            xp, yp = m(-8.6911,42.4061)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG MarinTG 
#            xp, yp = m(-9.074,39.586)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG NazareTG 
#            xp, yp = m(-9.367,39.354)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG PenicheTG 
#            xp, yp = m(-3.791,43.461)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG SantanderTG 
#            xp, yp = m(-8.888,37.949)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG SinesTG 
#            xp, yp = m(-1.67310,43.39840)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#17becf') # TG SocoaTG 
#            xp, yp = m(-8.726,42.243)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG VigoTG 
#            xp, yp = m(-8.77,42.601)
#            plt.scatter(xp, yp, s=100, alpha=1, c='#bcbd22') # TG VillagarciaTG 
            # Save and close 
            plt.savefig(plotname)
            plt.clf()

