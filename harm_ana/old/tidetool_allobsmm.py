#
# Script for DETIDING
# by Salish Sea MEOPAR
# https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis/raw/tip/compare_tides/Analysis8Components.ipynb
#
# imports
#%matplotlib inline
import matplotlib.pyplot as plt
#import matplotlib.pylab as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
#from salishsea_tools import tidetools
#from salishsea_tools import viz_tools
#from salishsea_tools import bathy_tools
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter 

path = '/work/ag15419/tmp/h_ana2/'
iniend_dates= '01/07/2017 - 31/12/2017'
dates_lab='20170701_20171231'
name = path
########################################################



# ============ ISPRA ===============

# STZ EX
#stations_mod = ['trieste','lampedusa','otranto']
#stations_obs = ['trieste','lampedusa','otranto']

# STZ EMODNET 2016 ALL (67)
#stations_obs = ['AjaccioTG-10min','AjaccioTG-1min','AjaccioTG-60min','AjaccioTG','AlcudiaTG','AlgecirasTG','AlmeriaTG','BarcelonaTG','CenturiTG-10min','CenturiTG-1min','CenturiTG-60min','CenturiTG','FormenteraTG','FosSurMerTG-10min','FosSurMerTG-1min','FosSurMerTG-60min','FosSurMerTG','GandiaTG','IbizaTG','IleRousseTG-10min','IleRousseTG-1min','IleRousseTG-60min','IleRousseTG','LaFigueiretteTG-10min','LaFigueiretteTG-1min','LaFigueiretteTG-60min','LaFigueiretteTG','MahonTG','MalagaTG','MarseilleTG-10min','MarseilleTG-60min','MarseilleTG','MelillaTG','MonacoTG-10min','MonacoTG-60min','MonacoTG','MotrilTG','NiceTG-10min','NiceTG-60min','NiceTG','PalmadeMallorcaTG','PortCamargueTG-5min','PortCamargueTG','PortFerreolTG-10min','PortFerreolTG-1min','PortFerreolTG-60min','PortFerreolTG','PortLaNouvelleTG-10min','PortLaNouvelleTG-1min','PortLaNouvelleTG-60min','PortLaNouvelleTG','PortVendresTG-10min','PortVendresTG-1min','PortVendresTG-60min','PortVendresTG','SaguntoTG','SeteTG-10min','SeteTG-1min','SeteTG-60min','SeteTG','SolenzaraTG-10min','SolenzaraTG-1min','SolenzaraTG-60min','SolenzaraTG','TarifaTG','ToulonTG','ValenciaTG']

# STZ EMODNET 2016 WORKING (60)
#stations_obs = ['AjaccioTG-10min','AjaccioTG-1min','AjaccioTG-60min','AjaccioTG','AlgecirasTG','AlmeriaTG','BarcelonaTG','CenturiTG-10min','CenturiTG-1min','CenturiTG-60min','CenturiTG','FormenteraTG','FosSurMerTG-10min','FosSurMerTG-1min','FosSurMerTG-60min','FosSurMerTG','GandiaTG','IbizaTG','IleRousseTG-10min','IleRousseTG-1min','IleRousseTG-60min','IleRousseTG','LaFigueiretteTG-10min','LaFigueiretteTG-1min','LaFigueiretteTG-60min','LaFigueiretteTG','MalagaTG','MarseilleTG-10min','MarseilleTG-60min','MarseilleTG','MelillaTG','MonacoTG-10min','MonacoTG-60min','MonacoTG','MotrilTG','PalmadeMallorcaTG','PortCamargueTG-5min','PortCamargueTG','PortFerreolTG-10min','PortFerreolTG-1min','PortFerreolTG-60min','PortFerreolTG','PortLaNouvelleTG-10min','PortLaNouvelleTG-1min','PortLaNouvelleTG-60min','PortLaNouvelleTG','PortVendresTG-10min','PortVendresTG-1min','PortVendresTG-60min','PortVendresTG','SaguntoTG','SeteTG-10min','SeteTG-1min','SeteTG-60min','SeteTG','SolenzaraTG-10min','SolenzaraTG-1min','SolenzaraTG-60min','SolenzaraTG','TarifaTG','ValenciaTG']
#stations_mod = ['AjaccioTG-10min','AjaccioTG-1min','AjaccioTG-60min','AjaccioTG','AlgecirasTG','AlmeriaTG','BarcelonaTG','CenturiTG-10min','CenturiTG-1min','CenturiTG-60min','CenturiTG','FormenteraTG','FosSurMerTG-10min','FosSurMerTG-1min','FosSurMerTG-60min','FosSurMerTG','GandiaTG','IbizaTG','IleRousseTG-10min','IleRousseTG-1min','IleRousseTG-60min','IleRousseTG','LaFigueiretteTG-10min','LaFigueiretteTG-1min','LaFigueiretteTG-60min','LaFigueiretteTG','MalagaTG','MarseilleTG-10min','MarseilleTG-60min','MarseilleTG','MelillaTG','MonacoTG-10min','MonacoTG-60min','MonacoTG','MotrilTG','PalmadeMallorcaTG','PortCamargueTG-5min','PortCamargueTG','PortFerreolTG-10min','PortFerreolTG-1min','PortFerreolTG-60min','PortFerreolTG','PortLaNouvelleTG-10min','PortLaNouvelleTG-1min','PortLaNouvelleTG-60min','PortLaNouvelleTG','PortVendresTG-10min','PortVendresTG-1min','PortVendresTG-60min','PortVendresTG','SaguntoTG','SeteTG-10min','SeteTG-1min','SeteTG-60min','SeteTG','SolenzaraTG-10min','SolenzaraTG-1min','SolenzaraTG-60min','SolenzaraTG','TarifaTG','ValenciaTG']

# STZ EMODNET 2016 HOURLY and RENAMED (25)
#stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortCamargue','PortFerreol','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
#stations_mod = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortCamargue','PortFerreol','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
#stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.Camargue','P.Ferreol','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']

# STZ ISPRA WORKING (16: 'ancona','carloforte','catania','civitavecchia','crotone','imperia','lampedusa','livorno','messina','napoli','ortona','otranto','palinuro','trieste','venezia','vieste')
#stations_obs = ['ancona','carloforte','catania','civitavecchia','crotone','imperia','lampedusa','livorno','messina','ortona','otranto','palinuro','trieste','venezia','vieste']
#stations_mod = ['ancona','carloforte','catania','civitavecchia','crotone','imperia','lampedusa','livorno','messina','ortona','otranto','palinuro','trieste','venezia','vieste']
#stations_mod2 = ['ancona','carloforte','catania','civitavecchia','crotone','imperia','lampedusa','livorno','messina','ortona','otranto','palinuro','trieste','venezia','vieste']
#stations_lab = ['Ancona','Carloforte','Catania','Civitavecchia','Crotone','Imperia','Lampedusa','Livorno','Messina','Ortona','Otranto','Palinuro','Trieste','Venezia','Vieste']
# x 2017
stations_obs = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
stations_mod = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
stations_mod2 = ['ancona','carloforte','catania','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
stations_lab = ['Ancona','Carloforte','Catania','Imperia','Lampedusa','Livorno','Messina','Ortona','Palinuro','Trieste','Venezia','Vieste']

# STZ ISPRA OK 
#stations_obs = ['ancona','carloforte','catania','crotone','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
#stations_mod = ['ancona','carloforte','catania','crotone','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
#stations_mod2 = ['ancona','carloforte','catania','crotone','imperia','lampedusa','livorno','messina','ortona','palinuro','trieste','venezia','vieste']
#stations_lab = ['Ancona','Carloforte','Catania','Crotone','Imperia','Lampedusa','Livorno','Messina','Ortona','Palinuro','Trieste','Venezia','Vieste']
#print("Stations list: ",stations)
##
numsta_obs=len(stations_obs)
numsta_mod=len(stations_mod)
numsta_mod2=len(stations_mod2)
#print("Stations num:",numsta )
##
###constants and fitting
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
#print("M2freq",M2freq,"radians per hour")
#print("S2freq",S2freq,"radians per hour")

### initial phase calculation

# M2
M2ft = 1.03669420605307
M2uvt = -12.5753366566839
#M2ft = 1.0
#M2uvt = 0.7238

# S2
S2ft = 1.000000
S2uvt = 0.000000

# N2
#N2ft = 1.037550
#N2uvt = 0.00464
N2ft = 1.0
N2uvt = 0.0

# K1
K1ft = 0.886345136935456
K1uvt = -0.4407793529972482
#K1ft = 1.0
#K1uvt = 0.0

# O1
O1ft = 0.813447400642605
O1uvt = -12.5067505117739
#O1ft = 1.0
#O1uvt = 0.0

# Q1
#Q1ft = 0.807594
#Q1uvt = -0.02993
Q1ft = 1.0
Q1uvt = 0.0

# K2
#K2ft = 0.7479278 
#K2uvt = 0.039694
K2ft = 1.0
K2uvt = 0.0

# P1
#P1ft = 1.000000
#P1uvt = 0.000000
P1ft = 1.0
P1uvt = 0.0

#print("Initial phases")
##

def octuple_obs(x, M2amp_obs, M2pha_obs, K1amp_obs, K1pha_obs, O1amp_obs, O1pha_obs, S2amp_obs, S2pha_obs, P1amp_obs, P1pha_obs, N2amp_obs, N2pha_obs, Q1amp_obs, Q1pha_obs, K2amp_obs, K2pha_obs):
    return (M2amp_obs*np.cos(M2freq*x-M2pha_obs*np.pi/180.)+
            K1amp_obs*np.cos(K1freq*x-K1pha_obs*np.pi/180.)+
            O1amp_obs*np.cos(O1freq*x-O1pha_obs*np.pi/180.)+
            S2amp_obs*np.cos(S2freq*x-S2pha_obs*np.pi/180.)+
            P1amp_obs*np.cos(P1freq*x-P1pha_obs*np.pi/180.)+
            N2amp_obs*np.cos(N2freq*x-N2pha_obs*np.pi/180.)+
            Q1amp_obs*np.cos(Q1freq*x-Q1pha_obs*np.pi/180.)+
            K2amp_obs*np.cos(K2freq*x-K2pha_obs*np.pi/180.))

def octuple_mod(x, M2amp_mod, M2pha_mod, K1amp_mod, K1pha_mod, O1amp_mod, O1pha_mod, S2amp_mod, S2pha_mod, P1amp_mod, P1pha_mod, N2amp_mod, N2pha_mod, Q1amp_mod, Q1pha_mod, K2amp_mod, K2pha_mod):
    return (M2amp_mod*np.cos(M2freq*x-M2pha_mod*np.pi/180.)+
            K1amp_mod*np.cos(K1freq*x-K1pha_mod*np.pi/180.)+
            O1amp_mod*np.cos(O1freq*x-O1pha_mod*np.pi/180.)+
            S2amp_mod*np.cos(S2freq*x-S2pha_mod*np.pi/180.)+
            P1amp_mod*np.cos(P1freq*x-P1pha_mod*np.pi/180.)+
            N2amp_mod*np.cos(N2freq*x-N2pha_mod*np.pi/180.)+
            Q1amp_mod*np.cos(Q1freq*x-Q1pha_mod*np.pi/180.)+
            K2amp_mod*np.cos(K2freq*x-K2pha_mod*np.pi/180.))

def octuple_mod2(x, M2amp_mod2, M2pha_mod2, K1amp_mod2, K1pha_mod2, O1amp_mod2, O1pha_mod2, S2amp_mod2, S2pha_mod2, P1amp_mod2, P1pha_mod2, N2amp_mod2, N2pha_mod2, Q1amp_mod2, Q1pha_mod2, K2amp_mod2, K2pha_mod2):
    return (M2amp_mod2*np.cos(M2freq*x-M2pha_mod2*np.pi/180.)+
            K1amp_mod2*np.cos(K1freq*x-K1pha_mod2*np.pi/180.)+
            O1amp_mod2*np.cos(O1freq*x-O1pha_mod2*np.pi/180.)+
            S2amp_mod2*np.cos(S2freq*x-S2pha_mod2*np.pi/180.)+
            P1amp_mod2*np.cos(P1freq*x-P1pha_mod2*np.pi/180.)+
            N2amp_mod2*np.cos(N2freq*x-N2pha_mod2*np.pi/180.)+
            Q1amp_mod2*np.cos(Q1freq*x-Q1pha_mod2*np.pi/180.)+
            K2amp_mod2*np.cos(K2freq*x-K2pha_mod2*np.pi/180.))


##

#allocate space for our arrays
M2_amp_obs=[]; M2_pha_obs=[]; K1_amp_obs=[]; K1_pha_obs=[]
O1_amp_obs=[]; O1_pha_obs=[]; S2_amp_obs=[]; S2_pha_obs=[]
P1_amp_obs=[]; P1_pha_obs=[]; N2_amp_obs=[]; N2_pha_obs=[]
Q1_amp_obs=[]; Q1_pha_obs=[]; K2_amp_obs=[]; K2_pha_obs=[]

M2_amp_mod=[]; M2_pha_mod=[]; K1_amp_mod=[]; K1_pha_mod=[]
O1_amp_mod=[]; O1_pha_mod=[]; S2_amp_mod=[]; S2_pha_mod=[]
P1_amp_mod=[]; P1_pha_mod=[]; N2_amp_mod=[]; N2_pha_mod=[]
Q1_amp_mod=[]; Q1_pha_mod=[]; K2_amp_mod=[]; K2_pha_mod=[]

M2_amp_mod2=[]; M2_pha_mod2=[]; K1_amp_mod2=[]; K1_pha_mod2=[]
O1_amp_mod2=[]; O1_pha_mod2=[]; S2_amp_mod2=[]; S2_pha_mod2=[]
P1_amp_mod2=[]; P1_pha_mod2=[]; N2_amp_mod2=[]; N2_pha_mod2=[]
Q1_amp_mod2=[]; Q1_pha_mod2=[]; K2_amp_mod2=[]; K2_pha_mod2=[]

#for stn in (0,1,2):
#for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60):
# EMODNET
#for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25):
# ISPRA
for stn in (0,1,2,3,4,5,6,7,8,9,10,11):
#for stn in (0,1,2,3,4,5,6,7):
    print('STZ:',stations_mod[stn])

    #fT1_obs = pd.read_csv(name+'obs_'+stations_obs[stn]+'.csv',sep=';',header=0,index_col='idNum', usecols=['idNum', 'station', 'year', 'month', 'day', 'hour', 'min', 'sec', 'value'])
    fT1_obs = pd.read_csv(name+'obs_'+stations_obs[stn]+'.csv',sep=';',usecols=['idNum', 'station', 'year', 'month', 'day', 'hour', 'minu', 'sec', 'value'])
    #years_obs=list(itemgetter(*fT1_obs.idNum)(fT1_obs.year)) 
    #year_obs=np.squeeze(fT1_obs.year[0,:])
    #print (year_obs)
    fT1_mod2 = NC.Dataset(name+stations_obs[stn]+'_mod_Tides4.nc','r')
    fT1_mod = NC.Dataset(name+stations_mod[stn]+'_mod_Tides8.nc','r')
    
    #print (name+stations_obs[stn]+'TG_mod_Tides4.nc')
    #print (name+stations_mod[stn]+'TG_mod_Tides8.nc')
    # ISPRA
    # prova = (datetime.datetime(2016, 7, 1, 0, 0, 0))
    # prova2 = (datetime.datetime(fT1_obs.year[:],fT1_obs.month[:],fT1_obs.day[:],fT1_obs.hour[:],fT1_obs.minu[:]))
    # prova3= (prova2-prova).total_seconds()
    #print (prova)
    #print (prova2)
    #print (prova3)
    #time_obs = (datetime.datetime(fT1_obs.year[:],fT1_obs.month[:],fT1_obs.day[:],fT1_obs.hour[:],fT1_obs.minu[:])-datetime.datetime(2016, 7, 1, 0, 0, 0)).total_seconds()
    #print ('time obs: ',time_obs)
    ssh_obs_arr= fT1_obs.value[:] 
    #print ('ssh: ',ssh_obs_arr)
    # NEMO:
    time_mod = (fT1_mod.variables["time_counter"][:]/3600.)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
    print('time mod:',time_mod)
    ssh_mod = fT1_mod.variables["sossheig"][:,0,0]*100.0 # want cm not meters
    #print(ssh_mod)
    # NEMO 2
    time_mod2 = (fT1_mod2.variables["time_counter"][:]/3600.)-582912.0 
    print ('time mod2:',time_mod2)
    ssh_mod2 = fT1_mod2.variables["sossheig"][:,0,0]*100.0
    #print(ssh_mod2)
    # times
    te_obs= ssh_obs_arr.shape[0]
    te_mod2 = ssh_mod2.shape[0]
    te_mod = ssh_mod.shape[0]
    #te_obs= 2232
    #
    ts_mod2=0
    ts_mod=0
    ts_obs=0

    # ISPRA
    time_obs=[]
    ssh_obs=[]
    #te_obs=4399
    for otime in range (0,te_obs):
        time_obs_el = (datetime.datetime(fT1_obs.year[otime],fT1_obs.month[otime],fT1_obs.day[otime],fT1_obs.hour[otime],fT1_obs.minu[otime])-datetime.datetime(2016, 7, 1, 0, 0, 0)).total_seconds()/3600 
        #time_obs.append(time_obs_el)
        # 
        num_of_nan=0
        if str(ssh_obs_arr[otime]) != 'nan' : # rm nan values from arrays
           #print ('ok:',ssh_obs_arr[otime])
           time_obs.append(time_obs_el)
           ssh_obs.append(float(ssh_obs_arr[otime]))
           #print ('obs:',time_obs[otime],ssh_obs[otime])
        else:
           num_of_nan=num_of_nan+1

    print ('# of nan = ',num_of_nan)
    #print ('time obs: ',time_obs)
    #print ('ssh obs:',ssh_obs)
    # Fit
    print("Fit mod2...")
    fitted_mod2, cov_mod2 = curve_fit(octuple_mod2,time_mod2[ts_mod2:te_mod2],ssh_mod2[ts_mod2:te_mod2])
    print ("Fit mod...")
    fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:te_mod],ssh_mod[ts_mod:te_mod])
    print ("Fit obs...")
    #valid = ~(np.isnan(ssh_obs)) 
    fitted_obs, cov_obs = curve_fit(octuple_obs,time_obs[ts_obs:te_obs],ssh_obs[ts_obs:te_obs])


    txtname=path+"pre_imm"+stations_mod[stn]+'_'+dates_lab+".txt"
    # Write arrays to file
    np.savetxt(txtname, (fitted_mod2,fitted_mod,fitted_obs), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")

     
    # Compute A and P in the choosen range
    if fitted_mod2[0] < 0:
        print('Am2_mod2<0')
        fitted_mod2[0] = -fitted_mod2[0]
        fitted_mod2[1] = fitted_mod2[1]+180
    M2_amp_mod2.append(fitted_mod2[0]*M2ft)
    pha_mod2 = fitted_mod2[1]+M2uvt

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    elif pha_mod2 < -360:
         pha_mod2=pha_mod2+360
         print('Pm2_mod2<-360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if pha_mod2 < 0:
        print('Pm2_mod2<0')
        pha_mod2 = pha_mod2+360

    fitted_mod2[1] = pha_mod2
    M2_pha_mod2.append(pha_mod2)
    
    if fitted_mod2[2] < 0:
        print (fitted_mod2[2])
        print('=> Ak1_mod2<0')
        fitted_mod2[2] = - fitted_mod2[2]
        fitted_mod2[3] = fitted_mod2[3] + 180
    K1_amp_mod2.append(fitted_mod2[2]*K1ft)
    pha_mod2 = fitted_mod2[3] + K1uvt
    if  pha_mod2 > 360:
        print('Pk1_mod2>360')
        pha_mod2 = pha_mod2-360
    if pha_mod2 < 0:
        print('Pk1_mod2<0')
        pha_mod2 = pha_mod2+360

    if pha_mod2 < 0:
        print('Pk1_mod2<0')
        pha_mod2 = pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    fitted_mod2[3] = pha_mod2
    K1_pha_mod2.append(pha_mod2)  
    
    if fitted_mod2[4] < 0:
        print('Ao1_mod2<0')
        fitted_mod2[4] = -fitted_mod2[4]
        fitted_mod2[5] = fitted_mod2[5]+180
    O1_amp_mod2.append(fitted_mod2[4]*O1ft)
    pha_mod2= fitted_mod2[5]+O1uvt
    if  pha_mod2 > 360:
        print('Po1_mod2>360')
        pha_mod2=pha_mod2-360
    elif pha_mod2 < 0:
        print('Po1_mod2<0')
        pha_mod2 = pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    fitted_mod2[5] = pha_mod2
    O1_pha_mod2.append(pha_mod2) 
    
    if fitted_mod2[6] < 0:
        print('As2_mod2<0')
        fitted_mod2[6] = -fitted_mod2[6]
        fitted_mod2[7] = fitted_mod2[7]+180
    S2_amp_mod2.append(fitted_mod2[6]*S2ft)
    pha_mod2= fitted_mod2[7]+S2uvt
    if  pha_mod2 > 360:
        print('Ps2_mod2>360')
        pha_mod2=pha_mod2-360
    elif pha_mod2 < -720:
        pha_mod2=pha_mod2+720
    elif pha_mod2 < -360:
        pha_mod2=pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if pha_mod2 < 0:
        print('Ps2_mod2<0')
        pha_mod2 = pha_mod2+360
    fitted_mod2[7] = pha_mod2
    S2_pha_mod2.append(pha_mod2) 
    
    
    #####
    if fitted_mod2[8] < 0:
            fitted_mod2[8] = -fitted_mod2[8]
            fitted_mod2[9] = fitted_mod2[9]+180
    P1_amp_mod2.append(fitted_mod2[8]*P1ft)
    pha_mod2= fitted_mod2[9]+P1uvt
    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
    elif pha_mod2 < 0:
        pha_mod2 = pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    fitted_mod2[9] = pha_mod2
    P1_pha_mod2.append(pha_mod2) 
    
    if fitted_mod2[10] < 0:
            fitted_mod2[10] = -fitted_mod2[10]
            fitted_mod2[11] = fitted_mod2[11]+180
    N2_amp_mod2.append(fitted_mod2[10]*N2ft)
    pha_mod2= fitted_mod2[11]+N2uvt
    if pha_mod2 > 720:
       pha_mod2 = pha_mod2-720
    elif pha_mod2 > 360:
        pha_mod2 = pha_mod2-360
    elif pha_mod2 < -360:
        pha_mod2 = pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')


    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    fitted_mod2[11] = pha_mod2
    N2_pha_mod2.append(pha_mod2) 
    
    if fitted_mod2[12] < 0:
            fitted_mod2[12] = -fitted_mod2[12]
            fitted_mod2[13] = fitted_mod2[13]+180
    Q1_amp_mod2.append(fitted_mod2[12]*Q1ft)
    pha_mod2= fitted_mod2[13]+Q1uvt
    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
    elif pha_mod2 < -360:
       pha_mod2=pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    fitted_mod2[13] = pha_mod2
    Q1_pha_mod2.append(pha_mod2) 
    
    if fitted_mod2[14] < 0:
            fitted_mod2[14] = -fitted_mod2[14]
            fitted_mod2[15] = fitted_mod2[15]+180
    K2_amp_mod2.append(fitted_mod2[14]*K2ft)
    pha_mod2= fitted_mod2[15]+K2uvt
    if pha_mod2 > 720:
       pha_mod2 = pha_mod2-720
    elif pha_mod2 > 360:
        pha_mod2 = pha_mod2-360
    elif pha_mod2 < -360:
        pha_mod2 = pha_mod2+360

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    fitted_mod2[15] = pha_mod2
    K2_pha_mod2.append(pha_mod2) 
    
    print('##########FIT OBS')
    print(fitted_mod2[:])
    #file.write('# FIT OBS')
    #np.savetxt('file_mod2mod.txt', zip(fitted_mod2,fitted_mod), fmt="%5.2f %5.2f")
    
    if fitted_mod[0] < 0:
        print('Am2_mod<0')
        fitted_mod[0] = -fitted_mod[0]
        fitted_mod[1] = fitted_mod[1]+180
    M2_amp_mod.append(fitted_mod[0]*M2ft)

    pha_mod = fitted_mod[1]+M2uvt
    if  pha_mod > 360:
        print('Pm2_mod>360')
        pha_mod=pha_mod-360

    if  pha_mod > 360:
        print('Pm2_mod>360')
        pha_mod=pha_mod-360

    elif pha_mod < -360:
        pha_mod=pha_mod+360

    if pha_mod < 0:
         print('Pm2_mod<0 ==> mod: ')
         pha_mod = pha_mod+360
         print (pha_mod)

    fitted_mod[1] = pha_mod
    M2_pha_mod.append(pha_mod)
    
    if fitted_mod[2] < 0:
        print('Ak1_mod<0')
        fitted_mod[2] = - fitted_mod[2]
        fitted_mod[3] = fitted_mod[3] + 180
    K1_amp_mod.append(fitted_mod[2]*K1ft)
    pha_mod = fitted_mod[3] + K1uvt
    if  pha_mod > 360:
        print('Pk1_mod>360')
        pha_mod = pha_mod-360

    if pha_mod < 0:
        print('Pk1_mod<0')
        pha_mod = pha_mod+360

    if pha_mod < 0:
        print('Pk1_mod<0')
        pha_mod = pha_mod+360

    fitted_mod[3] = pha_mod
    K1_pha_mod.append(pha_mod)
    
    if fitted_mod[4] < 0:
        print('Ao1_mod<0')
        fitted_mod[4] = -fitted_mod[4]
        fitted_mod[5] = fitted_mod[5]+180
    O1_amp_mod.append(fitted_mod[4]*O1ft)
    pha_mod= fitted_mod[5]+O1uvt
    if  pha_mod > 360:
        print('Po1_mod>360')
        pha_mod=pha_mod-360
    elif pha_mod < 0:
        print('Po1_mod<0')
        pha_mod = pha_mod+360

    fitted_mod[5] = pha_mod
    O1_pha_mod.append(pha_mod)
    
    if fitted_mod[6] < 0:
        print('As2_mod<0')
        fitted_mod[6] = -fitted_mod[6]
        fitted_mod[7] = fitted_mod[7]+180
    S2_amp_mod.append(fitted_mod[6]*S2ft)
    pha_mod= fitted_mod[7]+S2uvt

    if  pha_mod > 360:
        print('Ps2_mod>360')
        pha_mod=pha_mod-360
    elif pha_mod < -360:
        pha_mod=pha_mod+360

    if pha_mod < 0:
        print('Ps2_mod<0')
        pha_mod = pha_mod+360

    fitted_mod[7] = pha_mod
    S2_pha_mod.append(pha_mod)
    
    
    #####
    if fitted_mod[8] < 0:
            fitted_mod[8] = -fitted_mod[8]
            fitted_mod[9] = fitted_mod[9]+180
    P1_amp_mod.append(fitted_mod[8]*P1ft)
    pha_mod= fitted_mod[9]+P1uvt
    if  pha_mod > 360:
        pha_mod=pha_mod-360
    elif pha_mod < 0:
        pha_mod = pha_mod+360
    fitted_mod[9] = pha_mod
    P1_pha_mod.append(pha_mod)

    if fitted_mod[10] < 0:
            fitted_mod[10] = -fitted_mod[10]
            fitted_mod[11] = fitted_mod[11]+180
    N2_amp_mod.append(fitted_mod[10]*N2ft)
    pha_mod= fitted_mod[11]+N2uvt

    if pha_mod > 720:
       pha_mod = pha_mod-720
    elif pha_mod > 360:
        pha_mod = pha_mod-360
    elif pha_mod < -360:
        pha_mod = pha_mod+360

    if pha_mod > 360:
        pha_mod = pha_mod-360

    if pha_mod < 0:
       pha_mod = pha_mod+360

    fitted_mod[11] = pha_mod
    N2_pha_mod.append(pha_mod)
    
    if fitted_mod[12] < 0:
            fitted_mod[12] = -fitted_mod[12]
            fitted_mod[13] = fitted_mod[13]+180
    Q1_amp_mod.append(fitted_mod[12]*Q1ft)
    pha_mod= fitted_mod[13]+Q1uvt
    if  pha_mod > 360:
        pha_mod=pha_mod-360
    elif pha_mod < -360:
       pha_mod=pha_mod+360

    if pha_mod < 0:
        pha_mod = pha_mod+360
    fitted_mod[13] = pha_mod
    Q1_pha_mod.append(pha_mod)
    
    if fitted_mod[14] < 0:
            fitted_mod[14] = -fitted_mod[14]
            fitted_mod[15] = fitted_mod[15]+180
    K2_amp_mod.append(fitted_mod[14]*K2ft)
    pha_mod= fitted_mod[15]+K2uvt
    if pha_mod > 720:
       pha_mod = pha_mod-720
    elif pha_mod > 360:
        pha_mod = pha_mod-360
    elif pha_mod < -360:
        pha_mod = pha_mod+360

    if pha_mod < 0:
        pha_mod = pha_mod+360
    fitted_mod[15] = pha_mod
    K2_pha_mod.append(pha_mod)
    
    print('##########FIT MOD')
    print(fitted_mod[:])

##################
### OBS ###
##################
    if fitted_obs[0] < 0:
        print('Am2_obs<0')
        fitted_obs[0] = -fitted_obs[0]
        fitted_obs[1] = fitted_obs[1]+180
    M2_amp_obs.append(fitted_obs[0]*M2ft)

    pha_obs = fitted_obs[1]+M2uvt
    if  pha_obs > 360:
        print('Pm2_obs>360')
        pha_obs=pha_obs-360

    if  pha_obs > 360:
        print('Pm2_obs>360')
        pha_obs=pha_obs-360

    elif pha_obs < -360:
        pha_obs=pha_obs+360

    if pha_obs < 0:
         print('Pm2_obs<0 ==> mod: ')
         pha_obs = pha_obs+360
         print (pha_obs)

    fitted_obs[1] = pha_obs
    M2_pha_obs.append(pha_obs)

    if fitted_obs[2] < 0:
        print('Ak1_obs<0')
        fitted_obs[2] = - fitted_obs[2]
        fitted_obs[3] = fitted_obs[3] + 180
    K1_amp_obs.append(fitted_obs[2]*K1ft)
    pha_obs = fitted_obs[3] + K1uvt
    if  pha_obs > 360:
        print('Pk1_obs>360')
        pha_obs = pha_obs-360

    if pha_obs < 0:
        print('Pk1_obs<0')
        pha_obs = pha_obs+360

    if pha_obs < 0:
        print('Pk1_obs<0')
        pha_obs = pha_obs+360

    fitted_obs[3] = pha_obs
    K1_pha_obs.append(pha_obs)

    if fitted_obs[4] < 0:
        print('Ao1_obs<0')
        fitted_obs[4] = -fitted_obs[4]
        fitted_obs[5] = fitted_obs[5]+180
    O1_amp_obs.append(fitted_obs[4]*O1ft)
    pha_obs= fitted_obs[5]+O1uvt
    if  pha_obs > 360:
        print('Po1_obs>360')
        pha_obs=pha_obs-360
    elif pha_obs < 0:
        print('Po1_obs<0')
        pha_obs = pha_obs+360

    fitted_obs[5] = pha_obs
    O1_pha_obs.append(pha_obs)

    if fitted_obs[6] < 0:
        print('As2_obs<0')
        fitted_obs[6] = -fitted_obs[6]
        fitted_obs[7] = fitted_obs[7]+180
    S2_amp_obs.append(fitted_obs[6]*S2ft)
    pha_obs= fitted_obs[7]+S2uvt

    if  pha_obs > 360:
        print('Ps2_obs>360')
        pha_obs=pha_obs-360
    if pha_obs < -360:
        pha_obs=pha_obs+360
    if pha_obs < 0:
        print('Ps2_obs<0')
        pha_obs = pha_obs+360
    if pha_obs < 0:
        print('Ps2_obs<0')
        pha_obs = pha_obs+360
    fitted_obs[7] = pha_obs
    S2_pha_obs.append(pha_obs)


    #####
    if fitted_obs[8] < 0:
            fitted_obs[8] = -fitted_obs[8]
            fitted_obs[9] = fitted_obs[9]+180
    P1_amp_obs.append(fitted_obs[8]*P1ft)
    pha_obs= fitted_obs[9]+P1uvt
    if  pha_obs > 360:
        pha_obs=pha_obs-360
    elif pha_obs < 0:
        pha_obs = pha_obs+360
    fitted_obs[9] = pha_obs
    P1_pha_obs.append(pha_obs)

    if fitted_obs[10] < 0:
            fitted_obs[10] = -fitted_obs[10]
            fitted_obs[11] = fitted_obs[11]+180
    N2_amp_obs.append(fitted_obs[10]*N2ft)
    pha_obs= fitted_obs[11]+N2uvt

    if pha_obs > 720:
       pha_obs = pha_obs-720
    elif pha_obs > 360:
        pha_obs = pha_obs-360
    elif pha_obs < -360:
        pha_obs = pha_obs+360

    if pha_obs > 360:
        pha_obs = pha_obs-360

    if pha_obs < 0:
       pha_obs = pha_obs+360

    fitted_obs[11] = pha_obs
    N2_pha_obs.append(pha_obs)

    if fitted_obs[12] < 0:
            fitted_obs[12] = -fitted_obs[12]
            fitted_obs[13] = fitted_obs[13]+180
    Q1_amp_obs.append(fitted_obs[12]*Q1ft)
    pha_obs= fitted_obs[13]+Q1uvt
    if  pha_obs > 360:
        pha_obs=pha_obs-360
    elif pha_obs < -360:
       pha_obs=pha_obs+360

    if pha_obs < 0:
        pha_obs = pha_obs+360
    fitted_obs[13] = pha_obs
    Q1_pha_obs.append(pha_obs)

    if fitted_obs[14] < 0:
            fitted_obs[14] = -fitted_obs[14]
            fitted_obs[15] = fitted_obs[15]+180
    K2_amp_obs.append(fitted_obs[14]*K2ft)
    pha_obs= fitted_obs[15]+K2uvt
    if pha_obs > 720:
       pha_obs = pha_obs-720
    elif pha_obs > 360:
        pha_obs = pha_obs-360
    elif pha_obs < -360:
        pha_obs = pha_obs+360

    if pha_obs < 0:
        pha_obs = pha_obs+360
    fitted_obs[15] = pha_obs
    K2_pha_obs.append(pha_obs)

    print('##########FIT OBS')
    print(fitted_obs[:])


    txtname=path+"imm"+stations_mod[stn]+'_'+dates_lab+".txt"
    # Write arrays to file
    np.savetxt(txtname, (fitted_mod2,fitted_mod,fitted_obs), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")

print ("Plot...")

tidal_c=['M2','N2','K2','S2','K1','O1','Q1','P1']

for comp in ('M2','N2','K2','S2','K1','O1','Q1','P1'):

    #plt.figure(figsize=(24,10))
    ###
    #plt.subplot(2,1,1)
    nameA_obs=comp+'_amp_obs'
    nameA_mod=comp+'_amp_mod'
    nameA_mod2=comp+'_amp_mod2'

    nameP_obs=comp+'_pha_obs'
    nameP_mod=comp+'_pha_mod'
    nameP_mod2=comp+'_pha_mod2'

    # Diff computation
    #diff_A=globals()[nameA_mod]-globals()[nameA_obs]
    # diff Amp
    x_textA=[]
    y_textA=[]
    x_textA=np.array(globals()[nameA_obs])
    y_textA2=np.array(globals()[nameA_mod2])
    y_textA=np.array(globals()[nameA_mod])
    diffA_mo=y_textA-x_textA
    diffA_m2o=y_textA2-x_textA
    pdiffA_mo=(y_textA-x_textA)*100.0/x_textA
    pdiffA_m2o=(y_textA2-x_textA)*100.0/x_textA
    # diff Pha
    x_textP=[]
    y_textP=[]
    x_textP=np.array(globals()[nameP_obs])
    y_textP2=np.array(globals()[nameP_mod2])
    y_textP=np.array(globals()[nameP_mod])
    diffP_mo=y_textP-x_textP
    diffP_m2o=y_textP2-x_textP


    # Plot 1 ( A and P x stz )
    #plt.figure(figsize=(24,10))
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    # Amp 
    #plt.subplot(2,1,1)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(np.array(globals()[nameA_mod2]), 'g-o', label = 'Tides 4')
    plt.plot(stations_lab, np.array(globals()[nameA_mod]), 'r-o', label = 'Tides 8')
    plt.plot(np.array(globals()[nameA_obs]), '-bo', label = 'Obs')
    plt.title(comp+': Amplitude [cm] '+iniend_dates)
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.grid ()
    plt.ylabel ('Amplitude [cm]')
    # Amp diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(stations_lab,diffA_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(stations_lab, diffA_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.title(comp+': Model - Obs [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    plt.ylabel ('Model - Obs [cm]')
    # Add %
    #words=str(pdiffA_mo)
    #words=stations_lab
    i=0
    for word in pdiffA_mo:
         #plt.text(x_textA[i],diffA_m2o[i],word,fontsize=8,color = 'green')
         word_approx=round(word,1)
         string_to=str(word_approx)+'%'
         plt.text(stations_lab[i],diffA_mo[i]+.03,string_to,fontsize=10,color = 'black')
         i=i+1

    # Amp % diff
    #
    plt.savefig(path+comp+'_'+dates_lab+'_ADimm.jpg')
    ###
    # Pha
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    #plt.subplot(2,1,2)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(np.array(globals()[nameP_mod2]), 'g-o', label = 'Tides 4')
    plt.plot(stations_lab, np.array(globals()[nameP_mod]), 'r-o', label = 'Tides 8')
    plt.plot(np.array(globals()[nameP_obs]), '-bo', label = 'Obs')
    plt.title(comp+' Phase [deg] '+iniend_dates)
    plt.grid ()
    plt.ylim(0.0, 360.0)
    plt.legend( loc='upper left',fontsize = 'large' )
    # Pha diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(diffP_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(stations_lab,diffP_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.title(comp+' Model - Obs [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    #
    plt.savefig(path+comp+'_'+dates_lab+'_PDimm.jpg')

    # Polar Pha plot
    plt.figure(figsize=(12,12))
    plt.rc('font', size=10)
    #plt.title(comp+' Amplitudes (radius) [cm] and Phases (angle) [deg] ')
    ax=plt.subplot(111,projection='polar')
    #ax.set_rscale('log')
    #logA=log(radius)
    #logA1=log(radius1)
    #logA2=log(radius2)
    #radius=stations_lab
    theta=np.array(globals()[nameP_obs])
    theta1=np.array(globals()[nameP_mod])
    theta2=np.array(globals()[nameP_mod2])
    radius=np.array(globals()[nameA_obs])
    radius1=np.array(globals()[nameA_mod])
    radius2=np.array(globals()[nameA_mod2])
    plt.plot(theta,radius, 'bo', label = 'Obs')
    plt.plot(theta2,radius2, 'go', label = 'Tides 4')
    plt.plot(theta1,radius1, 'ro', label = 'Tides 8')
    #plt.plot(diffP_m2o,radius, 'go', label = 'Tides_4 - Obs')
    #plt.plot(diffP_mo,stations_lab, 'ro', label = 'Tides_8 - Obs')
    i=0
    for word in stations_lab:
        plt.text(theta[i],radius[i],word,fontsize=10,color = 'blue')
        plt.text(theta2[i],radius2[i],word,fontsize=10,color = 'green')
        plt.text(theta1[i],radius1[i],word,fontsize=10,color = 'red')
        i=i+1

    plt.legend( loc='upper left', fontsize = 'large' )
    plt.title(comp+' Amplitudes (radius) [cm] and Phases (angle) [deg] ')
    plt.savefig(path+comp+'_'+dates_lab+'_PolDimm.jpg')

    # Save ISPRA arrays
   
    IAmod2_name='I_'+nameA_mod2
    IAmod_name='I_'+nameA_mod
    IAobs_name='I_'+nameA_obs

    globals()[IAmod2_name]=np.array(globals()[nameA_mod2])
    globals()[IAmod_name]=np.array(globals()[nameA_mod])
    globals()[IAobs_name]=np.array(globals()[nameA_obs])

    I_stations_lab=stations_lab
    I_numsta=numsta_obs

    IPmod2_name='I_'+nameP_mod2
    IPmod_name='I_'+nameP_mod
    IPobs_name='I_'+nameP_obs


    globals()[IPmod2_name]=np.array(globals()[nameP_mod2])
    globals()[IPmod_name]=np.array(globals()[nameP_mod])
    globals()[IPobs_name]=np.array(globals()[nameP_obs])
  

# ======== EMODnet ========== 

# STZ EMODNET 2016 HOURLY and RENAMED (25)
#stations_mod2 = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortCamargue','PortFerreol','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
#stations_mod = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortCamargue','PortFerreol','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
#stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortCamargue','PortFerreol','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
#stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Formentera','FosSurMer','Gandia','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.Camargue','P.Ferreol','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']

# STZ EMODNET OK 2017 HOURLY and RENAMED 
stations_mod2 = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
stations_mod = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
stations_obs = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','PalmadeMallorca','PortLaNouvelle','PortVendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']
stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']

stations_lab = ['Ajaccio','Algeciras','Almeria','Barcelona','Centuri','Ibiza','IleRousse','LaFigueirette','Malaga','Marseille','Melilla','Monaco','Motril','P.deMallorca','P.LaNouvelle','P.Vendres','Sagunto','Sete','Solenzara','Tarifa','Valencia']

numsta_obs=len(stations_obs)
numsta_mod2=len(stations_mod2)
numsta_mod=len(stations_mod)
#print("Stations num:",numsta )
##
##
###constants and fitting
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
#print("M2freq",M2freq,"radians per hour")
#print("S2freq",S2freq,"radians per hour")

### initial phase calculation

# M2
M2ft = 1.03669420605307
M2uvt = -12.5753366566839
#M2ft = 1.0
#M2uvt = 0.7238

# S2
S2ft = 1.000000
S2uvt = 0.000000

# N2
#N2ft = 1.037550
#N2uvt = 0.00464
N2ft = 1.0
N2uvt = 0.0

# K1
K1ft = 0.886345136935456
K1uvt = -0.4407793529972482
#K1ft = 1.0
#K1uvt = 0.0

# O1
O1ft = 0.813447400642605
O1uvt = -12.5067505117739
#O1ft = 1.0
#O1uvt = 0.0

# Q1
#Q1ft = 0.807594
#Q1uvt = -0.02993
Q1ft = 1.0
Q1uvt = 0.0

# K2
#K2ft = 0.7479278 
#K2uvt = 0.039694
K2ft = 1.0
K2uvt = 0.0

# P1
#P1ft = 1.000000
#P1uvt = 0.000000
P1ft = 1.0
P1uvt = 0.0

#print("Initial phases")
##

def octuple_mod2(x, M2amp_mod2, M2pha_mod2, K1amp_mod2, K1pha_mod2, O1amp_mod2, O1pha_mod2, S2amp_mod2, S2pha_mod2, P1amp_mod2, P1pha_mod2, N2amp_mod2, N2pha_mod2, Q1amp_mod2, Q1pha_mod2, K2amp_mod2, K2pha_mod2):
    return (M2amp_mod2*np.cos(M2freq*x-M2pha_mod2*np.pi/180.)+
            K1amp_mod2*np.cos(K1freq*x-K1pha_mod2*np.pi/180.)+
            O1amp_mod2*np.cos(O1freq*x-O1pha_mod2*np.pi/180.)+
            S2amp_mod2*np.cos(S2freq*x-S2pha_mod2*np.pi/180.)+
            P1amp_mod2*np.cos(P1freq*x-P1pha_mod2*np.pi/180.)+
            N2amp_mod2*np.cos(N2freq*x-N2pha_mod2*np.pi/180.)+
            Q1amp_mod2*np.cos(Q1freq*x-Q1pha_mod2*np.pi/180.)+
            K2amp_mod2*np.cos(K2freq*x-K2pha_mod2*np.pi/180.))

def octuple_mod(x, M2amp_mod, M2pha_mod, K1amp_mod, K1pha_mod, O1amp_mod, O1pha_mod, S2amp_mod, S2pha_mod, P1amp_mod, P1pha_mod, N2amp_mod, N2pha_mod, Q1amp_mod, Q1pha_mod, K2amp_mod, K2pha_mod):
    return (M2amp_mod*np.cos(M2freq*x-M2pha_mod*np.pi/180.)+
            K1amp_mod*np.cos(K1freq*x-K1pha_mod*np.pi/180.)+
            O1amp_mod*np.cos(O1freq*x-O1pha_mod*np.pi/180.)+
            S2amp_mod*np.cos(S2freq*x-S2pha_mod*np.pi/180.)+
            P1amp_mod*np.cos(P1freq*x-P1pha_mod*np.pi/180.)+
            N2amp_mod*np.cos(N2freq*x-N2pha_mod*np.pi/180.)+
            Q1amp_mod*np.cos(Q1freq*x-Q1pha_mod*np.pi/180.)+
            K2amp_mod*np.cos(K2freq*x-K2pha_mod*np.pi/180.))

def octuple_obs(x, M2amp_obs, M2pha_obs, K1amp_obs, K1pha_obs, O1amp_obs, O1pha_obs, S2amp_obs, S2pha_obs, P1amp_obs, P1pha_obs, N2amp_obs, N2pha_obs, Q1amp_obs, Q1pha_obs, K2amp_obs, K2pha_obs):
    return (M2amp_obs*np.cos(M2freq*x-M2pha_obs*np.pi/180.)+
            K1amp_obs*np.cos(K1freq*x-K1pha_obs*np.pi/180.)+
            O1amp_obs*np.cos(O1freq*x-O1pha_obs*np.pi/180.)+
            S2amp_obs*np.cos(S2freq*x-S2pha_obs*np.pi/180.)+
            P1amp_obs*np.cos(P1freq*x-P1pha_obs*np.pi/180.)+
            N2amp_obs*np.cos(N2freq*x-N2pha_obs*np.pi/180.)+
            Q1amp_obs*np.cos(Q1freq*x-Q1pha_obs*np.pi/180.)+
            K2amp_obs*np.cos(K2freq*x-K2pha_obs*np.pi/180.))


##

#allocate space for our arrays
M2_amp_mod2=[]; M2_pha_mod2=[]; K1_amp_mod2=[]; K1_pha_mod2=[]
O1_amp_mod2=[]; O1_pha_mod2=[]; S2_amp_mod2=[]; S2_pha_mod2=[]
P1_amp_mod2=[]; P1_pha_mod2=[]; N2_amp_mod2=[]; N2_pha_mod2=[]
Q1_amp_mod2=[]; Q1_pha_mod2=[]; K2_amp_mod2=[]; K2_pha_mod2=[]

M2_amp_mod=[]; M2_pha_mod=[]; K1_amp_mod=[]; K1_pha_mod=[]
O1_amp_mod=[]; O1_pha_mod=[]; S2_amp_mod=[]; S2_pha_mod=[]
P1_amp_mod=[]; P1_pha_mod=[]; N2_amp_mod=[]; N2_pha_mod=[]
Q1_amp_mod=[]; Q1_pha_mod=[]; K2_amp_mod=[]; K2_pha_mod=[]

M2_amp_obs=[]; M2_pha_obs=[]; K1_amp_obs=[]; K1_pha_obs=[]
O1_amp_obs=[]; O1_pha_obs=[]; S2_amp_obs=[]; S2_pha_obs=[]
P1_amp_obs=[]; P1_pha_obs=[]; N2_amp_obs=[]; N2_pha_obs=[]
Q1_amp_obs=[]; Q1_pha_obs=[]; K2_amp_obs=[]; K2_pha_obs=[]

M2_Diff_A=[]
M2_Diff_P=[]
#for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60):
# EMODNET
#for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25):
for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20):
# ISPRA
#for stn in (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15):
    print(stations_mod[stn])

    fT1_mod2 = NC.Dataset(name+stations_mod2[stn]+'TG_mod_Tides4.nc','r')
    fT1_mod = NC.Dataset(name+stations_mod[stn]+'TG_mod_Tides8.nc','r')
    fT1_obs = NC.Dataset(name+stations_obs[stn]+'TG_obs.nc','r')
    #print (name+stations_mod2[stn]+'TG_mod_Tides4.nc')
    #print (name+stations_mod[stn]+'TG_mod_Tides8.nc')
    # mod:
    time_mod = (fT1_mod.variables["time_counter"][:]/3600.)-582912.0   # want hours (not seconds) from 1/1/2015 (not 1/1/1950)
    print('time mod:',time_mod)
    ssh_mod = fT1_mod.variables["sossheig"][:,0,0]*100.0 # want cm not meters
    #print(ssh_mod)
    # mod2
    time_mod2 = (fT1_mod2.variables["time_counter"][:]/3600.)-582912.0
    print ('time mod2:',time_mod2)
    ssh_mod2 = fT1_mod2.variables["sossheig"][:,0,0]*100.0
    #print(ssh_mod2)
    # TG obs
    time_obs = (fT1_obs.variables["TIME"][:]*24)-582912.0 # want hours (not days)
    print ('time obs:',time_obs)
    ssh_obs = fT1_obs.variables["sossheig"][:,0]*100.0
    print(ssh_obs)


    # times
    te_mod2 = ssh_mod2.shape[0]
    print (te_mod2)
    te_mod = ssh_mod.shape[0]
    print (te_mod)
    te_obs = ssh_obs.shape[0]
    print (te_obs)

    ts_mod2=0
    ts_mod=0
    ts_obs=0
    # Fit
    print("Fit mod2...")
    fitted_mod2, cov_mod2 = curve_fit(octuple_mod2,time_mod2[ts_mod2:te_mod2],ssh_mod2[ts_mod2:te_mod2])
    print ("Fit mod...")
    fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:te_mod],ssh_mod[ts_mod:te_mod])
    print ("Fit obs...")
    fitted_obs, cov_obs = curve_fit(octuple_obs,time_obs[ts_obs:te_obs],ssh_obs[ts_obs:te_obs])



    txtname=path+"pre_omm"+stations_mod[stn]+'_'+dates_lab+".txt"
    # Write arrays to file
    np.savetxt(txtname, (fitted_mod2,fitted_mod,fitted_obs), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")
    if fitted_mod2[0] < 0:
        print('Am2_mod2<0')
        fitted_mod2[0] = -fitted_mod2[0]
        fitted_mod2[1] = fitted_mod2[1]+180
    M2_amp_mod2.append(fitted_mod2[0]*M2ft)
    pha_mod2 = fitted_mod2[1]+M2uvt

    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
        print('Pm2_mod2>360')

    elif pha_mod2 < -360:
         pha_mod2=pha_mod2+360
         print('Pm2_mod2<-360')

    if pha_mod2 < 0:
        print('Pm2_mod2<0')
        pha_mod2 = pha_mod2+360

    fitted_mod2[1] = pha_mod2
    M2_pha_mod2.append(pha_mod2)

    if fitted_mod2[2] < 0:
        print (fitted_mod2[2])
        print('=> Ak1_mod2<0')
        fitted_mod2[2] = - fitted_mod2[2]
        fitted_mod2[3] = fitted_mod2[3] + 180
    K1_amp_mod2.append(fitted_mod2[2]*K1ft)
    pha_mod2 = fitted_mod2[3] + K1uvt
    if  pha_mod2 > 360:
        print('Pk1_mod2>360')
        pha_mod2 = pha_mod2-360
    elif pha_mod2 < 0:
        print('Pk1_mod2<0')
        pha_mod2 = pha_mod2+360
    fitted_mod2[3] = pha_mod2
    K1_pha_mod2.append(pha_mod2)

    if fitted_mod2[4] < 0:
        print('Ao1_mod2<0')
        fitted_mod2[4] = -fitted_mod2[4]
        fitted_mod2[5] = fitted_mod2[5]+180
    O1_amp_mod2.append(fitted_mod2[4]*O1ft)
    pha_mod2= fitted_mod2[5]+O1uvt
    if  pha_mod2 > 360:
        print('Po1_mod2>360')
        pha_mod2=pha_mod2-360
    elif pha_mod2 < 0:
        print('Po1_mod2<0')
        pha_mod2 = pha_mod2+360
    fitted_mod2[5] = pha_mod2
    O1_pha_mod2.append(pha_mod2)

    if fitted_mod2[6] < 0:
        print('As2_mod2<0')
        fitted_mod2[6] = -fitted_mod2[6]
        fitted_mod2[7] = fitted_mod2[7]+180
    S2_amp_mod2.append(fitted_mod2[6]*S2ft)
    pha_mod2= fitted_mod2[7]+S2uvt
    if  pha_mod2 > 360:
        print('ATTENZIONE: Ps2_mod2>360')
        pha_mod2=pha_mod2-360
    elif pha_mod2 < -720:
        pha_mod2=pha_mod2+720
    elif pha_mod2 < -360:
        pha_mod2=pha_mod2+360

    if  pha_mod2 > 360:
        print('ATTENZIONE: Ps2_mod2>360')
        pha_mod2=pha_mod2-360

    if pha_mod2 < 0:
        print('Ps2_mod2<0')
        pha_mod2 = pha_mod2+360
    fitted_mod2[7] = pha_mod2
    S2_pha_mod2.append(pha_mod2)


    #####
    # P1 tide 4
    if fitted_mod2[8] < 0:
            fitted_mod2[8] = -fitted_mod2[8]
            fitted_mod2[9] = fitted_mod2[9]+180
    P1_amp_mod2.append(fitted_mod2[8]*P1ft)
    pha_mod2= fitted_mod2[9]+P1uvt
    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
    if pha_mod2 < 0:
       pha_mod2 = pha_mod2+360
    if pha_mod2 < 0:
       pha_mod2 = pha_mod2+360
    if pha_mod2 < 0:
       pha_mod2 = pha_mod2+360
    fitted_mod2[9] = pha_mod2
    P1_pha_mod2.append(pha_mod2)

    # N2 tide 4
    if fitted_mod2[10] < 0:
            fitted_mod2[10] = -fitted_mod2[10]
            fitted_mod2[11] = fitted_mod2[11]+180
    N2_amp_mod2.append(fitted_mod2[10]*N2ft)
    pha_mod2= fitted_mod2[11]+N2uvt
    if pha_mod2 > 720:
       pha_mod2 = pha_mod2-720
    if pha_mod2 > 360:
       pha_mod2 = pha_mod2-360
    if pha_mod2 < -360:
       pha_mod2 = pha_mod2+360
    if pha_mod2 < 0:
       pha_mod2 = pha_mod2+360
    fitted_mod2[11] = pha_mod2
    N2_pha_mod2.append(pha_mod2)

    if fitted_mod2[12] < 0:
            fitted_mod2[12] = -fitted_mod2[12]
            fitted_mod2[13] = fitted_mod2[13]+180
    Q1_amp_mod2.append(fitted_mod2[12]*Q1ft)
    pha_mod2= fitted_mod2[13]+Q1uvt
    if  pha_mod2 > 360:
        pha_mod2=pha_mod2-360
    elif pha_mod2 < -360:
       pha_mod2=pha_mod2+360

    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    fitted_mod2[13] = pha_mod2
    Q1_pha_mod2.append(pha_mod2)

    # K2 tide 4
    if fitted_mod2[14] < 0:
            fitted_mod2[14] = -fitted_mod2[14]
            fitted_mod2[15] = fitted_mod2[15]+180
    K2_amp_mod2.append(fitted_mod2[14]*K2ft)
    pha_mod2= fitted_mod2[15]+K2uvt
    if pha_mod2 > 720:
       pha_mod2 = pha_mod2-720
    elif pha_mod2 > 360:
        pha_mod2 = pha_mod2-360
    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    if pha_mod2 < 0:
        pha_mod2 = pha_mod2+360
    fitted_mod2[15] = pha_mod2
    K2_pha_mod2.append(pha_mod2)

    print('##########FIT MOD2')
    print(fitted_mod2[:])
    #file.write('# FIT OBS')
    #np.savetxt('file_mod2mod.txt', zip(fitted_mod2,fitted_mod), fmt="%5.2f %5.2f")
    if fitted_mod[0] < 0:
        print('Am2_mod<0')
        fitted_mod[0] = -fitted_mod[0]
        fitted_mod[1] = fitted_mod[1]+180
    M2_amp_mod.append(fitted_mod[0]*M2ft)

    pha_mod = fitted_mod[1]+M2uvt
    if  pha_mod > 360:
        print('Pm2_mod>360')
        pha_mod=pha_mod-360

    elif pha_mod < -360:
        pha_mod=pha_mod+360

    if pha_mod < 0:
         print('Pm2_mod<0 ==> mod: ')
         pha_mod = pha_mod+360
         print (pha_mod)

    fitted_mod[1] = pha_mod
    M2_pha_mod.append(pha_mod)

    if fitted_mod[2] < 0:
        print('Ak1_mod<0')
        fitted_mod[2] = - fitted_mod[2]
        fitted_mod[3] = fitted_mod[3] + 180
    K1_amp_mod.append(fitted_mod[2]*K1ft)
    pha_mod = fitted_mod[3] + K1uvt
    if  pha_mod > 360:
        print('Pk1_mod>360')
        pha_mod = pha_mod-360
    elif pha_mod < 0:
        print('Pk1_mod<0')
        pha_mod = pha_mod+360
    fitted_mod[3] = pha_mod
    K1_pha_mod.append(pha_mod)

    # O1 tide 8
    if fitted_mod[4] < 0:
        print('Ao1_mod<0')
        fitted_mod[4] = -fitted_mod[4]
        fitted_mod[5] = fitted_mod[5]+180
    O1_amp_mod.append(fitted_mod[4]*O1ft)
    pha_mod= fitted_mod[5]+O1uvt
    if  pha_mod > 360:
        print('Po1_mod>360')
        pha_mod=pha_mod-360
    elif pha_mod < 0:
        print('Po1_mod<0')
        pha_mod = pha_mod+360
    if  pha_mod > 360:
        print('Po1_mod>360')
        pha_mod=pha_mod-360
    if  pha_mod > 360:
        print('Po1_mod>360')
        pha_mod=pha_mod-360

    fitted_mod[5] = pha_mod
    O1_pha_mod.append(pha_mod)

    if fitted_mod[6] < 0:
        print('As2_mod<0')
        fitted_mod[6] = -fitted_mod[6]
        fitted_mod[7] = fitted_mod[7]+180
    S2_amp_mod.append(fitted_mod[6]*S2ft)
    pha_mod= fitted_mod[7]+S2uvt

    if  pha_mod > 360:
        print('Ps2_mod>360')
        pha_mod=pha_mod-360
    elif pha_mod < -360:
        pha_mod=pha_mod+360

    if  pha_mod > 360:
        print('Ps2_mod>360')
        pha_mod=pha_mod-360

    if pha_mod < 0:
        print('Ps2_mod<0')
        pha_mod = pha_mod+360

    fitted_mod[7] = pha_mod
    S2_pha_mod.append(pha_mod)


    #####
    # P1 tide 8
    if fitted_mod[8] < 0:
            fitted_mod[8] = -fitted_mod[8]
            fitted_mod[9] = fitted_mod[9]+180
    P1_amp_mod.append(fitted_mod[8]*P1ft)
    pha_mod= fitted_mod[9]+P1uvt
    if  pha_mod > 360:
        pha_mod=pha_mod-360
    if pha_mod < 0:
       pha_mod = pha_mod+360
    fitted_mod[9] = pha_mod
    P1_pha_mod.append(pha_mod)

    if fitted_mod[10] < 0:
            fitted_mod[10] = -fitted_mod[10]
            fitted_mod[11] = fitted_mod[11]+180
    N2_amp_mod.append(fitted_mod[10]*N2ft)
    pha_mod= fitted_mod[11]+N2uvt
    if pha_mod > 720:
       pha_mod = pha_mod-720
    elif pha_mod > 360:
        pha_mod = pha_mod-360
    elif pha_mod < -360:
        pha_mod = pha_mod+360

    if pha_mod < 0:
        pha_mod = pha_mod+360
    fitted_mod[11] = pha_mod
    N2_pha_mod.append(pha_mod)

    # Q1 tide 8
    if fitted_mod[12] < 0:
            fitted_mod[12] = -fitted_mod[12]
            fitted_mod[13] = fitted_mod[13]+180
    Q1_amp_mod.append(fitted_mod[12]*Q1ft)
    pha_mod= fitted_mod[13]+Q1uvt
    if  pha_mod > 360:
        pha_mod=pha_mod-360
    elif pha_mod < -360:
       pha_mod=pha_mod+360
    if  pha_mod > 360:
        pha_mod=pha_mod-360

    if pha_mod < 0:
        pha_mod = pha_mod+360
    fitted_mod[13] = pha_mod
    Q1_pha_mod.append(pha_mod)

    # K2 tide 8
    if fitted_mod[14] < 0:
            fitted_mod[14] = -fitted_mod[14]
            fitted_mod[15] = fitted_mod[15]+180
    K2_amp_mod.append(fitted_mod[14]*K2ft)
    pha_mod= fitted_mod[15]+K2uvt
    if pha_mod > 720:
       pha_mod = pha_mod-720
    if pha_mod > 360:
       pha_mod = pha_mod-360
    if pha_mod < 0:
       pha_mod = pha_mod+360
    if pha_mod < 0:
        pha_mod = pha_mod+360
    if pha_mod < 0:
        pha_mod = pha_mod+360

    fitted_mod[15] = pha_mod
    K2_pha_mod.append(pha_mod)

    print('##########FIT MOD')
    print(fitted_mod[:])


    # M2 obs
    if fitted_obs[0] < 0:
        print('Am2_obs<0')
        fitted_obs[0] = -fitted_obs[0]
        fitted_obs[1] = fitted_obs[1]+180
    M2_amp_obs.append(fitted_obs[0]*M2ft)

    pha_obs = fitted_obs[1]+M2uvt
    if  pha_obs > 360:
        print('Pm2_obs>360')
        pha_obs=pha_obs-360
    if pha_obs < -360:
       pha_obs=pha_obs+360
    if pha_obs < 0:
         print('Pm2_obs<0 ==> mod: ')
         pha_obs = pha_obs+360
         print (pha_obs)
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360

    fitted_obs[1] = pha_obs
    M2_pha_obs.append(pha_obs)

    # K1 obs
    if fitted_obs[2] < 0:
        print('Ak1_obs<0')
        fitted_obs[2] = - fitted_obs[2]
        fitted_obs[3] = fitted_obs[3] + 180
    K1_amp_obs.append(fitted_obs[2]*K1ft)
    pha_obs = fitted_obs[3] + K1uvt
    if  pha_obs > 360:
        print('Pk1_obs>360')
        pha_obs = pha_obs-360
    elif pha_obs < 0:
        print('Pk1_obs<0')
        pha_obs = pha_obs+360
    fitted_obs[3] = pha_obs
    K1_pha_obs.append(pha_obs)

    # O1 obs
    if fitted_obs[4] < 0:
        print('Ao1_obs<0')
        fitted_obs[4] = -fitted_obs[4]
        fitted_obs[5] = fitted_obs[5]+180
    O1_amp_obs.append(fitted_obs[4]*O1ft)
    pha_obs= fitted_obs[5]+O1uvt
    if  pha_obs > 720:
        pha_obs=pha_obs-720
    if  pha_obs > 360:
        print('Po1_obs>360')
        pha_obs=pha_obs-360
    if  pha_obs < 0:
        print('Po1_obs<0')
        pha_obs = pha_obs+360
    if  pha_obs > 360:
        print('Po1_obs>360')
        pha_obs=pha_obs-360
    if pha_obs < 0:
        print('Po1_obs<0')
        pha_obs = pha_obs+360
    if pha_obs < 0:
        print('Po1_obs<0')
        pha_obs = pha_obs+360

    fitted_obs[5] = pha_obs
    O1_pha_obs.append(pha_obs)

    # S2 obs
    if fitted_obs[6] < 0:
        print('As2_obs<0')
        fitted_obs[6] = -fitted_obs[6]
        fitted_obs[7] = fitted_obs[7]+180
    S2_amp_obs.append(fitted_obs[6]*S2ft)
    pha_obs= fitted_obs[7]+S2uvt

    if  pha_obs > 360:
        print('ATTENZIONE: Ps2_obs>360')
        pha_obs=pha_obs-360
    elif pha_obs < -360:
        pha_obs=pha_obs+360

    if  pha_obs > 360:
        print('ATTENZIONE: Ps2_obs>360')
        pha_obs=pha_obs-360

    if  pha_obs > 360:
        print('ATTENZIONE: Ps2_obs>360')
        pha_obs=pha_obs-360

    if pha_obs < 0:
        print('Ps2_obs<0')
        pha_obs = pha_obs+360

    fitted_obs[7] = pha_obs
    S2_pha_obs.append(pha_obs)


    #####

    # P1 obs
    if fitted_obs[8] < 0:
            fitted_obs[8] = -fitted_obs[8]
            fitted_obs[9] = fitted_obs[9]+180
    P1_amp_obs.append(fitted_obs[8]*P1ft)
    pha_obs= fitted_obs[9]+P1uvt
    if  pha_obs > 360:
        pha_obs=pha_obs-360
    elif pha_obs < 0:
        pha_obs = pha_obs+360
    fitted_obs[9] = pha_obs
    P1_pha_obs.append(pha_obs)

    # N2 obs
    if fitted_obs[10] < 0:
            fitted_obs[10] = -fitted_obs[10]
            fitted_obs[11] = fitted_obs[11]+180
    N2_amp_obs.append(fitted_obs[10]*N2ft)
    pha_obs= fitted_obs[11]+N2uvt
    if pha_obs > 720:
       pha_obs = pha_obs-720
    elif pha_obs > 360:
         pha_obs = pha_obs-360

    if pha_obs < -360:
       pha_obs = pha_obs+360

    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360

    fitted_obs[11] = pha_obs
    N2_pha_obs.append(pha_obs)


    # Q1 obs
    if fitted_obs[12] < 0:
            fitted_obs[12] = -fitted_obs[12]
            fitted_obs[13] = fitted_obs[13]+180
    Q1_amp_obs.append(fitted_obs[12]*Q1ft)
    pha_obs= fitted_obs[13]+Q1uvt
    if  pha_obs > 360:
        pha_obs=pha_obs-360
    elif pha_obs < -360:
       pha_obs=pha_obs+360

    if pha_obs < 0:
        pha_obs = pha_obs+360
    fitted_obs[13] = pha_obs
    Q1_pha_obs.append(pha_obs)

    # K2 obs
    if fitted_obs[14] < 0:
            fitted_obs[14] = -fitted_obs[14]
            fitted_obs[15] = fitted_obs[15]+180
    K2_amp_obs.append(fitted_obs[14]*K2ft)
    pha_obs= fitted_obs[15]+K2uvt
    if pha_obs > 720:
       pha_obs = pha_obs-720
    elif pha_obs > 360:
        pha_obs = pha_obs-360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    if pha_obs < 0:
       pha_obs = pha_obs+360
    fitted_obs[15] = pha_obs
    K2_pha_obs.append(pha_obs)

    print('##########FIT OBS')
    print(fitted_obs[:])
    #M2_Diff_A[stn]=fitted_obs[:]-fitted_mod[:]



    txtname=path+"omm_"+stations_mod[stn]+'_'+dates_lab+"_mod.txt"
    # Write arrays to file
    #mod2_arr=M2_amp_mod2+M2_pha_mod2+K1_amp_mod2+K1_pha_mod2+O1_amp_mod2+O1_pha_mod2+S2_amp_mod2+S2_pha_mod2+P1_amp_mod2+P1_pha_mod2+N2_amp_mod2+N2_pha_mod2+Q1_amp_mod2+Q1_pha_mod2+K2_amp_mod2+K2_pha_mod2
    #mod_arr=M2_amp_mod+M2_pha_mod+K1_amp_mod+K1_pha_mod+O1_amp_mod+O1_pha_mod+S2_amp_mod+S2_pha_mod+P1_amp_mod+P1_pha_mod+N2_amp_mod+N2_pha_mod+Q1_amp_mod+Q1_pha_mod+K2_amp_mod+K2_pha_mod
    #print(mod2_arr)
    #mod2_arr.shape
    #np.squeeze(mod_arr).shape
    ##print(mod2_arr)
    np.savetxt(txtname, (fitted_mod2,fitted_mod,fitted_obs), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")


tidal_c=['M2','N2','K2','S2','K1','O1','Q1','P1']


for comp in ('M2','N2','K2','S2','K1','O1','Q1','P1'):
    #plt.figure(figsize=(24,10))
    ###
    #plt.subplot(2,1,1)
    nameA_obs=comp+'_amp_obs'
    nameA_mod=comp+'_amp_mod'
    nameA_mod2=comp+'_amp_mod2'

    nameP_obs=comp+'_pha_obs'
    nameP_mod=comp+'_pha_mod'
    nameP_mod2=comp+'_pha_mod2'

    # Diff computation
    #diff_A=globals()[nameA_mod]-globals()[nameA_obs]
    # diff Amp
    x_textA=[]
    y_textA=[]
    diffA_mo=[]
    diffA_m2o=[]
    pdiffA_mo=[]
    pdiffA_m2o=[]
    x_textA=np.array(globals()[nameA_obs])
    y_textA2=np.array(globals()[nameA_mod2])
    y_textA=np.array(globals()[nameA_mod])
    diffA_mo=y_textA-x_textA
    diffA_m2o=y_textA2-x_textA
    pdiffA_mo=(y_textA-x_textA)*100.0/x_textA
    pdiffA_m2o=(y_textA2-x_textA)*100.0/x_textA
    # diff Pha
    x_textP=[]
    y_textP=[]
    x_textP=np.array(globals()[nameP_obs])
    y_textP2=np.array(globals()[nameP_mod2])
    y_textP=np.array(globals()[nameP_mod])
    diffP_mo=y_textP-x_textP
    diffP_m2o=y_textP2-x_textP


    # Plot 1 ( A and P x stz )
    #plt.figure(figsize=(24,10))
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    # Amp 
    #plt.subplot(2,1,1)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(np.array(globals()[nameA_mod2]), 'g-o', label = 'Tides 4')
    plt.plot(stations_lab, np.array(globals()[nameA_mod]), 'r-o', label = 'Tides 8')
    plt.plot(np.array(globals()[nameA_obs]), '-bo', label = 'Obs')
    plt.title(comp+': Amplitude [cm] '+iniend_dates)
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.grid ()
    plt.ylabel ('Amplitude [cm]')
    # Amp diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(stations_lab,diffA_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(stations_lab, diffA_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.title(comp+': Model - Obs [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    plt.ylabel ('Model - Obs [cm]')
    # Add %
    #words=str(pdiffA_mo)
    #words=stations_lab
    i=0
    for word in pdiffA_mo:
         #plt.text(x_textA[i],diffA_m2o[i],word,fontsize=8,color = 'green')
         word_approx=round(word,1)
         string_to=str(word_approx)+'%'
         plt.text(stations_lab[i],diffA_mo[i]+.03,string_to,fontsize=10,color = 'black')
         i=i+1
#    for i in enumerate(stations_lab):
#        plt.text(stations_lab[i],diffA_mo[i]+.03,words[i],fontsize=8,color = 'red')

    plt.savefig(path+comp+'_'+dates_lab+'_ADomm.jpg')


    # Amp % diff
#    plt.subplot(3,1,3)
#    plt.plot(stations_lab,pdiffA_m2o, 'g-o', label = '(Tides_4 - Obs)*100/Obs')
#    plt.plot(stations_lab, pdiffA_mo, 'r-o', label = '(Tides_8 - Obs)*100/Obs')
#    plt.legend( loc='upper left' )
#    plt.title(comp+': (Model - Obs)*100/Obs [%]')
#    plt.grid ()
#    plt.axhline(y=0)
#    plt.ylabel ('(Model - Obs)*100/Obs [%]')
   ###
    # Pha
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    #plt.subplot(2,1,2)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(np.array(globals()[nameP_mod2]), 'g-o', label = 'Tides 4')
    plt.plot(stations_lab, np.array(globals()[nameP_mod]), 'r-o', label = 'Tides 8')
    plt.plot(np.array(globals()[nameP_obs]), '-bo', label = 'Obs')
    plt.title(comp+' Phase [deg] '+iniend_dates)
    plt.grid ()
    plt.ylabel ('Phase [deg]')
    plt.ylim(0.0, 360.0)
    plt.legend( loc='upper left',fontsize = 'large')
    # Pha diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(diffP_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(stations_lab,diffP_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.title(comp+' Model - Obs [deg] ')
    plt.grid ()
    plt.ylabel ('Model - Obs [deg]')
    plt.axhline(y=0, color='black')
    #
    plt.savefig(path+comp+'_'+dates_lab+'_PDomm.jpg')

    # Polar Pha plot
    plt.figure(figsize=(12,12))
    plt.rc('font', size=10)
    #plt.ylim(0.0, 5.0)
    ax=plt.subplot(111,projection='polar')
    #ax.set_rscale('log')
    #logA=log(radius)
    #logA1=log(radius1)
    #logA2=log(radius2)
    #radius=stations_lab
    #ax.set_rorigin(0)
    #plt.ylim(0.0, 5.0)
    theta=np.array(globals()[nameP_obs])
    theta1=np.array(globals()[nameP_mod])
    theta2=np.array(globals()[nameP_mod2])
    radius=np.array(globals()[nameA_obs])
    radius1=np.array(globals()[nameA_mod])
    radius2=np.array(globals()[nameA_mod2])
    plt.plot(theta,radius, 'bo', label = 'Obs')
    plt.plot(theta2,radius2, 'go', label = 'Tides 4')
    plt.plot(theta1,radius1, 'ro', label = 'Tides 8')
    #plt.plot(diffP_m2o,radius, 'go', label = 'Tides_4 - Obs')
    #plt.plot(diffP_mo,stations_lab, 'ro', label = 'Tides_8 - Obs')
    i=0
    for word in stations_lab:
        plt.text(theta[i],radius[i],word,fontsize=10,color = 'blue')
        plt.text(theta2[i],radius2[i],word,fontsize=10,color = 'green')
        plt.text(theta1[i],radius1[i],word,fontsize=10,color = 'red')
        i=i+1
    plt.title(comp+' Amplitudes (radius) [cm] and Phases (angle) [deg] ')
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.savefig(path+comp+'_'+dates_lab+'_PolDomm.jpg')

# ============ ISPRA + EMODnet


for comp in ('M2','N2','K2','S2','K1','O1','Q1','P1'):

  # Plot only T8 for N2, K2, P1, Q1
  if comp != 'M2' and comp != 'S2' and comp != 'K1' and comp != 'O1' : 

    nameA_obs=comp+'_amp_obs'
    nameA_mod=comp+'_amp_mod'

    nameP_obs=comp+'_pha_obs'
    nameP_mod=comp+'_pha_mod'
 
    I_nameA_obs='I_'+comp+'_amp_obs'
    I_nameA_mod='I_'+comp+'_amp_mod'

    I_nameP_obs='I_'+comp+'_pha_obs'
    I_nameP_mod='I_'+comp+'_pha_mod'

    TOT_A_obs=np.append(globals()[nameA_obs],globals()[I_nameA_obs])
    TOT_A_obs=np.append(np.array(globals()[nameA_obs]),np.array(globals()[I_nameA_obs]))    
    TOT_A_mod=np.append(np.array(globals()[nameA_mod]),np.array(globals()[I_nameA_mod]))    

    TOT_P_obs=np.append(np.array(globals()[nameP_obs]),np.array(globals()[I_nameP_obs]))    
    TOT_P_mod=np.append(np.array(globals()[nameP_mod]),np.array(globals()[I_nameP_mod]))

    
    TOT_stations_lab=np.append(stations_lab,I_stations_lab)
    lab=0
    for stz_nm in TOT_stations_lab:
        TOT_stations_lab[lab]=stz_nm[:3]
        lab=lab+1
        #print (stz_nm)

    #print (TOT_stations_lab)


    # Diff and diff_mean computation
    #diff_A=globals()[TOT_A_mod]-globals()[TOT_A_obs]
    # diff Amp
    x_textA=[]
    y_textA=[]
    x_textA=TOT_A_obs
    y_textA=TOT_A_mod
    diffA_mo=y_textA-x_textA
    pdiffA_mo=(y_textA-x_textA)*100.0/x_textA
    mean_diffA=np.mean(diffA_mo)

    # diff Pha
    x_textP=[]
    y_textP=[]
    x_textP=TOT_P_obs
    y_textP=TOT_P_mod
    diffP_mo=y_textP-x_textP
    mean_diffP=np.mean(diffP_mo)

    # Plot 1 ( A and P x stz )
    #plt.figure(figsize=(24,10))
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    # Amp 
    #plt.subplot(2,1,1)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(TOT_stations_lab, np.array(TOT_A_mod), 'r-o', label = 'Tides 8')
    plt.plot(TOT_A_obs, '-bo', label = 'Obs')
    plt.title(comp+': Amplitude [cm] '+iniend_dates)
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.grid ()
    plt.ylabel ('Amplitude [cm]')
    # Amp diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(TOT_stations_lab, diffA_mo, 'r-o', label = 'Tides_8 - Obs')
    #plt.legend( loc='upper left',fontsize = 'large' )
    plt.title(comp+': Model - Obs [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    App_mean_diffA=round(mean_diffA,2)
    plt.axhline(y=mean_diffA, color='red',label = '<Amp Diff>='+str(App_mean_diffA)+'[cm]')
    plt.ylabel ('Model - Obs [cm]')
    plt.legend( loc='upper left',fontsize = 'large' )
    # Add %
    #words=str(pdiffA_mo)
    #words=TOT_stations_lab
    i=0
    for word in pdiffA_mo:
         #plt.text(x_textA[i],diffA_m2o[i],word,fontsize=8,color = 'green')
         word_approx=round(word,1)
         string_to=str(word_approx)+'%'
         plt.text(TOT_stations_lab[i],diffA_mo[i]+.03,string_to,fontsize=10,color = 'black')
         i=i+1
#    for i in enumerate(TOT_stations_lab):
#        plt.text(TOT_stations_lab[i],diffA_mo[i]+.03,words[i],fontsize=8,color = 'red')

    plt.savefig(path+comp+'_'+dates_lab+'_Aallmm.jpg')
  
   ###
    # Pha
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    #plt.subplot(2,1,2)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(TOT_stations_lab, TOT_P_mod, 'r-o', label = 'Tides 8')
    plt.plot(TOT_P_obs, '-bo', label = 'Obs')
    plt.title(comp+' Phase [deg] '+iniend_dates)
    plt.grid ()
    plt.ylabel ('Phase [deg]')
    plt.ylim(0.0, 360.0)
    plt.legend( loc='upper left',fontsize = 'large')
    # Pha diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(TOT_stations_lab,diffP_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.title(comp+' Model - Obs [deg] ')
    plt.grid ()
    plt.ylabel ('Model - Obs [deg]')
    App_mean_diffP=round(mean_diffP,2)
    plt.axhline(y=mean_diffP, color='red',label = '<Pha Diff>='+str(App_mean_diffP)+'[deg]')
    plt.axhline(y=0, color='black')
    plt.legend( loc='upper left',fontsize = 'large' )
    #
    plt.savefig(path+comp+'_'+dates_lab+'_Pallmm.jpg')


######## LIN REG


   # Plot 2 ( Lin reg mod Vs obs x A and P x stz )

    plt.figure(figsize=(6,12))
    plt.rc('font', size=12)
    #
    plt.subplot(2,1,1)
    #plt.plot(np.array(TOT_A_obs),np.array(TOT_A_mod2),'go',label = 'Tides 4')
    #plt.plot(np.array(TOT_A_obs), np.array(TOT_A_mod),'ro',label = 'Tides 8')
    plt.title(comp+' Amplitude [cm] '+'('+iniend_dates+')')
    #plt.legend( loc='upper left' )
    plt.grid ()
    #bottom_x, top_x = plt.xlim()
    #bottom_y, top_y = plt.ylim()
    #top=max(top_x,top_y)
    #plt.ylim(0.0, top+1)
    #plt.xlim(0.0, top+1)
    #plt.plot([0.0,100], [0.0,100], 'k-', color = 'black')

    # Arrays defn
    x_text=[]
    y_text=[]
    y_text2=[]
    x_text=np.array(TOT_A_obs)
    #print ('x text: ',x_text)
    #y_text2=np.array(TOT_A_mod2)
    #print ('y text2: ',y_text2)
    y_text=np.array(TOT_A_mod)
    #print ('y text: ',y_text)  


    top=np.maximum(x_text,y_text)
    top=max(top[:])

    print('top: ',top)
    # Linear regression
    # Mod T8
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_text,y_text)
    m_A=[]
    q_A=[]
    fitted_A=[]
    cov_A=[]
    def line_A(x, m_A, q_A):
        return (m_A*x+q_A)
    fitted_A, cov_A = curve_fit(line_A,x_text[:],y_text[:])
    m_A=fitted_A[0]
    q_A=fitted_A[1]
    rx=np.linspace(0.0,top)
    retta_A=m_A*rx+q_A
    m_A_approx=round(m_A,2)
    #lr_leg_str='( slope='+str(m_A_approx)+')'
    #print ('y text: ',y_text)
    print (cov_A)
    #perr = np.abs(np.sqrt(np.diag(cov_A)))
    perr = np.abs(np.diag(cov_A))
    print(perr)
    m_Ae_approx=round(perr[0],2)
    r_A=round(r_value,2)
    #lr_leg_str='( slope='+str(m_A_approx)+'+-'+str(m_Ae_approx)+')'
    lr_leg_str='( slope='+str(m_A_approx)+'; R2='+str(r_A)+')'

    # Arrays defn
    x_text=[]
    #y_text=[]
    #y_text2=[]
    x_text=np.array(TOT_A_obs)
    #print ('x text: ',x_text)
    #y_text2=np.array(TOT_A_mod2)
    #print ('y text2: ',y_text2)
    #y_text=np.array(TOT_A_mod)
    #print ('y text: ',y_text)

    # Mod T4
    #slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x_text,y_text2)
    #m_A2=[]
    #q_A2=[]
    #fitted_A2=[]
    #cov_a2=[]
    #def line_A2(x, m_A2, q_A2):
    #    return (m_A2*x+q_A2)
    #fitted_A2, cov_A2 = curve_fit(line_A2,x_text[:],y_text2[:])
    #m_A2=fitted_A2[0]
    #q_A2=fitted_A2[1]
    #rx2=np.linspace(0.0,top)
    #retta_A2=m_A2*rx2+q_A2
    #m_A2_approx=round(m_A2,2)
    #lr_leg_str2='( slope='+str(m_A2_approx)+')'
    #print ('y_text2: ',y_text2)
    #print (cov_A2)
    #perr2 = np.abs(np.sqrt(np.diag(cov_A2)))
    #perr2 = np.abs(np.diag(cov_A2))
    #print(perr2)
    #m_A2e_approx=round(perr2[0],2)
    #r_A2=round(r_value2,2)
    #lr_leg_str2='( slope='+str(m_A2_approx)+'+-'+str(m_A2e_approx)+')'
    #lr_leg_str2='( slope='+str(m_A2_approx)+'; R2='+str(r_A2)+')'


    # Plot lines
    #plt.plot(rx2,retta_A2,color = 'green')
    plt.plot(rx,retta_A,color = 'red')
    plt.plot([0.0,top], [0.0,top], 'k-', color = 'black')
    #plt.plot(np.array(TOT_A_obs),np.array(TOT_A_mod2),'go',label = 'Tides 4 '+lr_leg_str2)
    plt.plot(np.array(TOT_A_obs), np.array(TOT_A_mod),'ro',label = 'Tides 8 '+lr_leg_str)

    plt.xlabel ('OBS Amplitude [cm]')
    plt.ylabel ('MOD Amplitude [cm]')
    plt.legend( loc='upper left' )

    # point label (stn names)
    i=0
    for word in TOT_stations_lab:
        #plt.text(x_text[i]+.03,y_text2[i]+.03,word,fontsize=8,color = 'green')
        plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=8,color = 'red')
        i=i+1

    ### Pha linear reg
    plt.subplot(2,1,2)
    #plt.plot(np.array(TOT_P_obs),np.array(TOT_P_mod2),'go', label = 'Tides 4')
    #plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 ')
    plt.title(comp+' Phase [deg] '+'('+iniend_dates+')')
    plt.grid ()
    plt.ylim(0.0, 360.0)
    plt.xlim(0.0, 360.0)
    plt.plot([0.0, 360.0], [0.0, 360.0], 'k-', color = 'black')
    x_text=[]
    y_text=[]
    #y_text2=[]
    x_text=np.array(TOT_P_obs)
    #y_text2=np.array(TOT_P_mod2)
    y_text=np.array(TOT_P_mod)
    #slopeP2, interceptP2, r_valueP2, p_valueP2, std_errP2 = stats.linregress(x_text,y_text2)
    slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
    #rx2=np.linspace(0.0,360.0)
    rx=np.linspace(0.0,360.0)
    #retta_P2=slopeP2*rx2+interceptP2
    retta_P=slopeP*rx+interceptP
    #plt.plot(rx2,retta_P2,color = 'green')
    plt.plot(rx,retta_P,color = 'red')
    #lr_leg_str2='( slope='+str(round(slopeP2,2))+'; R2='+str(round(r_valueP2,2))+')'
    #lr_leg_str='( slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
    #lr_leg_str2='( slope='+str(round(slopeP2,2))+'; R2='+str(round(r_valueP2,2))+')'
    lr_leg_str='( slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
    #plt.plot(np.array(TOT_P_obs),np.array(TOT_P_mod2),'go', label = 'Tides 4'+lr_leg_str2)
    plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 '+lr_leg_str)

    plt.xlabel ('OBS Phase [deg]')
    plt.ylabel ('MOD Phase [deg]')
    plt.legend( loc='upper left' )
    i=0
    for word in TOT_stations_lab:
        #plt.text(x_text[i]+.03,y_text2[i]+.03,word,fontsize=8,color = 'green')
        plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=8,color = 'red')
        i=i+1

    plt.savefig(path+comp+'_'+dates_lab+'_all_lr.jpg')


    # Compute RMSmisfit
    tot_stn_number=len(TOT_stations_lab)
    d_foreman=np.sqrt((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod*np.cos(TOT_P_mod)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod*np.sin(TOT_P_mod))**2))
    rms_misfit=np.sqrt((1/(2*tot_stn_number))*(np.sum((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod*np.cos(TOT_P_mod)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod*np.sin(TOT_P_mod))**2))))
    print (comp+' RMSmisfit: '+str(rms_misfit))
    print (comp+' d_foreman: '+str(d_foreman))

  else: 
  # T8 and T4 constituents
    nameA_obs=comp+'_amp_obs'
    nameA_mod=comp+'_amp_mod'
    nameA_mod2=comp+'_amp_mod2'

    nameP_obs=comp+'_pha_obs'
    nameP_mod=comp+'_pha_mod'
    nameP_mod2=comp+'_pha_mod2'

    I_nameA_obs='I_'+comp+'_amp_obs'
    I_nameA_mod='I_'+comp+'_amp_mod'
    I_nameA_mod2='I_'+comp+'_amp_mod2'

    I_nameP_obs='I_'+comp+'_pha_obs'
    I_nameP_mod='I_'+comp+'_pha_mod'
    I_nameP_mod2='I_'+comp+'_pha_mod2'

    TOT_A_obs=np.append(globals()[nameA_obs],globals()[I_nameA_obs])
    TOT_A_obs=np.append(np.array(globals()[nameA_obs]),np.array(globals()[I_nameA_obs]))
    TOT_A_mod=np.append(np.array(globals()[nameA_mod]),np.array(globals()[I_nameA_mod]))
    TOT_A_mod2=np.append(np.array(globals()[nameA_mod2]),np.array(globals()[I_nameA_mod2]))

    TOT_P_obs=np.append(np.array(globals()[nameP_obs]),np.array(globals()[I_nameP_obs]))
    TOT_P_mod=np.append(np.array(globals()[nameP_mod]),np.array(globals()[I_nameP_mod]))
    TOT_P_mod2=np.append(np.array(globals()[nameP_mod2]),np.array(globals()[I_nameP_mod2]))


    TOT_stations_lab=np.append(stations_lab,I_stations_lab)
    lab=0
    for stz_nm in TOT_stations_lab:
        TOT_stations_lab[lab]=stz_nm[:3]
        lab=lab+1
        #print (stz_nm)

    #print (TOT_stations_lab)
    # Diff and mean diff computation
    #diff_A=globals()[TOT_A_mod]-globals()[TOT_A_obs]
    # diff Amp
    x_textA=[]
    y_textA=[]
    x_textA=TOT_A_obs
    y_textA2=TOT_A_mod2
    y_textA=TOT_A_mod
    diffA_mo=y_textA-x_textA
    diffA_m2o=y_textA2-x_textA
    pdiffA_mo=(y_textA-x_textA)*100.0/x_textA
    pdiffA_m2o=(y_textA2-x_textA)*100.0/x_textA
    mean_diffA2=np.mean(diffA_m2o)
    mean_diffA=np.mean(diffA_mo)
    

    # diff Pha
    x_textP=[]
    y_textP=[]
    x_textP=TOT_P_obs
    y_textP2=TOT_P_mod2
    y_textP=TOT_P_mod
    diffP_mo=y_textP-x_textP
    diffP_m2o=y_textP2-x_textP
    mean_diffP2=np.mean(diffP_m2o)
    mean_diffP=np.mean(diffP_mo)


    # Plot 1 ( A and P x stz )
    #plt.figure(figsize=(24,10))
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    # Amp 
    #plt.subplot(2,1,1)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(TOT_A_mod2, 'g-o', label = 'Tides 4')
    plt.plot(TOT_stations_lab, np.array(TOT_A_mod), 'r-o', label = 'Tides 8')
    plt.plot(TOT_A_obs, '-bo', label = 'Obs')
    plt.title(comp+': Amplitude [cm] '+iniend_dates)
    plt.legend( loc='upper left',fontsize = 'large' )
    plt.grid ()
    plt.ylabel ('Amplitude [cm]')
    # Amp diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(TOT_stations_lab,diffA_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(TOT_stations_lab, diffA_mo, 'r-o', label = 'Tides_8 - Obs')
    #plt.legend( loc='upper left',fontsize = 'large' )
    plt.title(comp+': Model - Obs [cm] ')
    plt.grid ()
    plt.axhline(y=0, color='black')
    plt.ylabel ('Model - Obs [cm]')
    App_mean_diffA2=round(mean_diffA2,2)
    App_mean_diffA=round(mean_diffA,2)
    plt.axhline(y=mean_diffA2, color='green', label = '<Amp Diff T4>='+str(App_mean_diffA2)+'[cm]')
    plt.axhline(y=mean_diffA, color='red',label = '<Amp Diff T8>='+str(App_mean_diffA)+'[cm]')
    plt.legend( loc='upper left',fontsize = 'large' )
    # Add %
    #words=str(pdiffA_mo)
    #words=TOT_stations_lab
    i=0
    for word in pdiffA_mo:
         #plt.text(x_textA[i],diffA_m2o[i],word,fontsize=8,color = 'green')
         word_approx=round(word,1)
         string_to=str(word_approx)+'%'
         plt.text(TOT_stations_lab[i],diffA_mo[i]+.03,string_to,fontsize=10,color = 'black')
         i=i+1
#    for i in enumerate(TOT_stations_lab):
#        plt.text(TOT_stations_lab[i],diffA_mo[i]+.03,words[i],fontsize=8,color = 'red')

    plt.savefig(path+comp+'_'+dates_lab+'_Aallmm.jpg')

   ###
    # Pha
    plt.figure(figsize=(24,10))
    plt.rc('font', size=10)
    #plt.subplot(2,1,2)
    plt.subplot(2,1,1)
    plt.xticks(fontsize=10)
    plt.plot(TOT_P_mod2, 'g-o', label = 'Tides 4')
    plt.plot(TOT_stations_lab, TOT_P_mod, 'r-o', label = 'Tides 8')
    plt.plot(TOT_P_obs, '-bo', label = 'Obs')
    plt.title(comp+' Phase [deg] '+iniend_dates)
    plt.grid ()
    plt.ylabel ('Phase [deg]')
    plt.ylim(0.0, 360.0)
    plt.legend( loc='upper left',fontsize = 'large')
    # Pha diff
    plt.subplot(2,1,2)
    plt.xticks(fontsize=10)
    plt.plot(diffP_m2o, 'g-o', label = 'Tides_4 - Obs')
    plt.plot(TOT_stations_lab,diffP_mo, 'r-o', label = 'Tides_8 - Obs')
    plt.title(comp+' Model - Obs [deg] ')
    plt.grid ()
    plt.ylabel ('Model - Obs [deg]')
    plt.axhline(y=0, color='black')
    App_mean_diffP2=round(mean_diffP2,2)
    App_mean_diffP=round(mean_diffP,2)
    plt.axhline(y=mean_diffP2, color='red', label = '<Pha Diff T4>='+str(App_mean_diffP2)+'[deg]')
    plt.axhline(y=mean_diffP, color='green', label = '<Pha Diff T8>='+str(App_mean_diffP)+'[deg]')
    plt.legend( loc='upper left',fontsize = 'large' )
    #
    plt.savefig(path+comp+'_'+dates_lab+'_Pallmm.jpg')


######## LIN REG


   # Plot 2 ( Lin reg mod Vs obs x A and P x stz )

    plt.figure(figsize=(6,12))
    plt.rc('font', size=12)
    #
    plt.subplot(2,1,1)
    #plt.plot(np.array(TOT_A_obs),np.array(TOT_A_mod2),'go',label = 'Tides 4')
    #plt.plot(np.array(TOT_A_obs), np.array(TOT_A_mod),'ro',label = 'Tides 8')
    plt.title(comp+' Amplitude [cm] '+'('+iniend_dates+')')
    #plt.legend( loc='upper left' )
    plt.grid ()
    #bottom_x, top_x = plt.xlim()
    #bottom_y, top_y = plt.ylim()
    #top=max(top_x,top_y)
    #plt.ylim(0.0, top+1)
    #plt.xlim(0.0, top+1)
    #plt.plot([0.0,100], [0.0,100], 'k-', color = 'black')

    # Arrays defn
    x_text=[]
    y_text=[]
    y_text2=[]
    x_text=np.array(TOT_A_obs)
    #print ('x text: ',x_text)
    y_text2=np.array(TOT_A_mod2)
    #print ('y text2: ',y_text2)
    y_text=np.array(TOT_A_mod)
    #print ('y text: ',y_text)  


    top=np.maximum(x_text,y_text,y_text2)
    top=max(top[:])

    print('top: ',top)
    # Linear regression
    # Mod T8
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_text,y_text)
    m_A=[]
    q_A=[]
    fitted_A=[]
    cov_A=[]
    def line_A(x, m_A, q_A):
        return (m_A*x+q_A)
    fitted_A, cov_A = curve_fit(line_A,x_text[:],y_text[:])
    m_A=fitted_A[0]
    q_A=fitted_A[1]
    rx=np.linspace(0.0,top)
    retta_A=m_A*rx+q_A
    m_A_approx=round(m_A,2)
    #lr_leg_str='( slope='+str(m_A_approx)+')'
    #print ('y text: ',y_text)
    print (cov_A)
    #perr = np.abs(np.sqrt(np.diag(cov_A)))
    perr = np.abs(np.diag(cov_A))
    print(perr)
    m_Ae_approx=round(perr[0],2)
    r_A=round(r_value,2)
    #lr_leg_str='( slope='+str(m_A_approx)+'+-'+str(m_Ae_approx)+')'
    lr_leg_str='( slope='+str(m_A_approx)+'; R2='+str(r_A)+')'

    # Arrays defn
    x_text=[]
    #y_text=[]
    y_text2=[]
    x_text=np.array(TOT_A_obs)
    #print ('x text: ',x_text)
    y_text2=np.array(TOT_A_mod2)
    #print ('y text2: ',y_text2)
    #y_text=np.array(TOT_A_mod)
    #print ('y text: ',y_text)

    # Mod T4
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x_text,y_text2)
    m_A2=[]
    q_A2=[]
    fitted_A2=[]
    cov_a2=[]
    def line_A2(x, m_A2, q_A2):
        return (m_A2*x+q_A2)
    fitted_A2, cov_A2 = curve_fit(line_A2,x_text[:],y_text2[:])
    m_A2=fitted_A2[0]
    q_A2=fitted_A2[1]
    rx2=np.linspace(0.0,top)
    retta_A2=m_A2*rx2+q_A2
    m_A2_approx=round(m_A2,2)
    #lr_leg_str2='( slope='+str(m_A2_approx)+')'
    #print ('y_text2: ',y_text2)
    print (cov_A2)
    #perr2 = np.abs(np.sqrt(np.diag(cov_A2)))
    perr2 = np.abs(np.diag(cov_A2))
    print(perr2)
    m_A2e_approx=round(perr2[0],2)
    r_A2=round(r_value2,2)
    #lr_leg_str2='( slope='+str(m_A2_approx)+'+-'+str(m_A2e_approx)+')'
    lr_leg_str2='( slope='+str(m_A2_approx)+'; R2='+str(r_A2)+')'


    # Plot lines
    plt.plot(rx2,retta_A2,color = 'green')
    plt.plot(rx,retta_A,color = 'red')
    plt.plot([0.0,top], [0.0,top], 'k-', color = 'black')
    plt.plot(np.array(TOT_A_obs),np.array(TOT_A_mod2),'go',label = 'Tides 4 '+lr_leg_str2)
    plt.plot(np.array(TOT_A_obs), np.array(TOT_A_mod),'ro',label = 'Tides 8 '+lr_leg_str)

    plt.xlabel ('OBS Amplitude [cm]')
    plt.ylabel ('MOD Amplitude [cm]')
    plt.legend( loc='upper left' )

    # point label (stn names)
    i=0
    for word in TOT_stations_lab:
        plt.text(x_text[i]+.03,y_text2[i]+.03,word,fontsize=8,color = 'green')
        plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=8,color = 'red')
        i=i+1

    ### Pha linear reg
    plt.subplot(2,1,2)
    #plt.plot(np.array(TOT_P_obs),np.array(TOT_P_mod2),'go', label = 'Tides 4')
    #plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 ')
    plt.title(comp+' Phase [deg] '+'('+iniend_dates+')')
    plt.grid ()
    plt.ylim(0.0, 360.0)
    plt.xlim(0.0, 360.0)
    plt.plot([0.0, 360.0], [0.0, 360.0], 'k-', color = 'black')
    x_text=[]
    y_text=[]
    y_text2=[]
    x_text=np.array(TOT_P_obs)
    y_text2=np.array(TOT_P_mod2)
    y_text=np.array(TOT_P_mod)
    slopeP2, interceptP2, r_valueP2, p_valueP2, std_errP2 = stats.linregress(x_text,y_text2)
    slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
    rx2=np.linspace(0.0,360.0)
    rx=np.linspace(0.0,360.0)
    retta_P2=slopeP2*rx2+interceptP2
    retta_P=slopeP*rx+interceptP
    plt.plot(rx2,retta_P2,color = 'green')
    plt.plot(rx,retta_P,color = 'red')
    #lr_leg_str2='( slope='+str(round(slopeP2,2))+'; R2='+str(round(r_valueP2,2))+')'
    #lr_leg_str='( slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
    lr_leg_str2='( slope='+str(round(slopeP2,2))+'; R2='+str(round(r_valueP2,2))+')'
    lr_leg_str='( slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
    plt.plot(np.array(TOT_P_obs),np.array(TOT_P_mod2),'go', label = 'Tides 4'+lr_leg_str2)
    plt.plot(np.array(TOT_P_obs), np.array(TOT_P_mod),'ro',label = 'Tides 8 '+lr_leg_str)

    plt.xlabel ('OBS Phase [deg]')
    plt.ylabel ('MOD Phase [deg]')
    plt.legend( loc='upper left' )
    i=0
    for word in TOT_stations_lab:
        plt.text(x_text[i]+.03,y_text2[i]+.03,word,fontsize=8,color = 'green')
        plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=8,color = 'red')
        i=i+1

    plt.savefig(path+comp+'_'+dates_lab+'_all_lr.jpg')

    # Compute RMSmisfit
    tot_stn_number=len(TOT_stations_lab)
    d_foreman=np.sqrt((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod*np.cos(TOT_P_mod)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod*np.sin(TOT_P_mod))**2))
    rms_misfit=np.sqrt((1/(2*tot_stn_number))*(np.sum((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod*np.cos(TOT_P_mod)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod*np.sin(TOT_P_mod))**2))))
    print (comp+' RMSmisfit: '+str(rms_misfit))
    print (comp+' d_foreman: '+str(d_foreman))

    d_foreman2=np.sqrt((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod2*np.cos(TOT_P_mod2)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod2*np.sin(TOT_P_mod2))**2))
    rms_misfit2=np.sqrt((1/(2*tot_stn_number))*(np.sum((TOT_A_obs*np.cos(TOT_P_obs)-(TOT_A_mod2*np.cos(TOT_P_mod2)))**2-(TOT_A_obs*np.sin(TOT_P_obs)-(TOT_A_mod2*np.sin(TOT_P_mod2))**2))))
    print (comp+' RMSmisfit2: '+str(rms_misfit2))
    print (comp+' d_foreman2: '+str(d_foreman2))



