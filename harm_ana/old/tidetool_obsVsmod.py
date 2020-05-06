#
# Script for DETIDING
# by Salish Sea MEOPAR
# https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis/raw/tip/compare_tides/Analysis8Components.ipynb
#
#print("Setting the environment...")
# imports
#%matplotlib inline
import matplotlib.pyplot as plt
#import matplotlib.pylab as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
#from salishsea_tools import tidetools
#from salishsea_tools import viz_tools
#from salishsea_tools import bathy_tools
import collections
import pandas as pd
import csv
import math

#from __future__ import division
#print("..Done!")
###
#### pathname for data - all of the tide runs are stored in this directory
####path = '/data/nsoontie/MEOPAR/SalishSea/results/tides/'
path = '/work/ag15419/tmp/h_ana/'
#print("Work dir is",path )
###
####the run we want to analyze
####runname = 'corr15'
runname = 'simu_tides8'
#print("Runname is",runname )
###
####joining the two string together
#name = path +runname +'/'
name = path
###
#print("Name:",name)
###
#stations_mod = ['lampedusa','otranto','trieste']
#stations_obs = ['lampedusa','otranto','trieste']
stations_obs = ['MotrilTG']
stations_mod = ['MotrilTG']
#print("Stations list: ",stations)
###
numsta_obs=len(stations_obs)
numsta_mod=len(stations_mod)
#print("Stations num:",numsta )
###
####constants and fitting
#### M2
M2freq = 28.984106 # degrees per hour
M2freq = M2freq*np.pi/180. # radians per hour
####K1
K1freq = 15.041069*np.pi/180.
####O1
O1freq = 13.943036*np.pi/180.
####S2
S2freq = 30.000002*np.pi/180.
####P1
P1freq = 14.958932*np.pi/180.
####N2
N2freq = 28.439730*np.pi/180.
####Q1
Q1freq = 13.398661*np.pi/180.
####K2
K2freq = 30.082138*np.pi/180.
###
#print("M2freq",M2freq,"radians per hour")
#print("S2freq",S2freq,"radians per hour")

#### initial phase calculation
#M2ft = 1.036581
#M2uvt = 0.009395
M2ft = 1.0
M2uvt = 0.0


S2ft = 1.000000
S2uvt = 0.000000

#N2ft = 1.037550
#N2uvt = 0.00464
N2ft = 1.0
N2uvt = 0.0

#K1ft = 0.882850
#K1uvt = 0.022005
K1ft = 1.0
K1uvt = 0.0

#O1ft = 0.807594
#O1uvt = -0.02993
O1ft = 1.0
O1uvt = 0.0

Q1ft = 0.807594
Q1uvt = -0.02993

K2ft = 0.7479278 
K2uvt = 0.039694

P1ft = 1.000000
P1uvt = 0.000000

#print("Initial phases")
#print("M2 ft and uvt:",M2ft,M2uvt)
#print("S2 ft and uvt:",S2ft,S2uvt)
###
#def double(x, M2amp, M2pha, S2amp, S2pha):
#    return (M2amp*np.cos(M2freq*x-M2pha*np.pi/180.)+
#            S2amp*np.cos(S2freq*x-S2pha*np.pi/180.))

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



###
#fig, ax = plt.subplots(1,1,figsize=(12,5))

#allocate space for our arrays
M2_amp_obs=[]; M2_pha_obs=[]; K1_amp_obs=[]; K1_pha_obs=[]
O1_amp_obs=[]; O1_pha_obs=[]; S2_amp_obs=[]; S2_pha_obs=[]
P1_amp_obs=[]; P1_pha_obs=[]; N2_amp_obs=[]; N2_pha_obs=[]
Q1_amp_obs=[]; Q1_pha_obs=[]; K2_amp_obs=[]; K2_pha_obs=[]

M2_amp_mod=[]; M2_pha_mod=[]; K1_amp_mod=[]; K1_pha_mod=[]
O1_amp_mod=[]; O1_pha_mod=[]; S2_amp_mod=[]; S2_pha_mod=[]
P1_amp_mod=[]; P1_pha_mod=[]; N2_amp_mod=[]; N2_pha_mod=[]
Q1_amp_mod=[]; Q1_pha_mod=[]; K2_amp_mod=[]; K2_pha_mod=[]


#file=open("amp_pha.txt","w")

for stn in (0,0):
    print(stations_mod[stn])
    #file.write("# "+stations_mod[stn])
    #print(stations)
    fT1_obs = NC.Dataset(name+stations_obs[stn]+'_obs.nc','r')
    fT1_mod = NC.Dataset(name+stations_mod[stn]+'_mod_Tides8.nc','r')
    print (name+stations_obs[stn]+'_obs.nc')
    print (name+stations_mod[stn]+'_mod.nc')
    # NEMO:
    time_mod = fT1_mod.variables["time_counter"][:]/3600. # want hours not seconds
    print('time mod:',time_mod)
    ssh_mod = fT1_mod.variables["sossheig"][:,0,0]
    print(ssh_mod)
    # TG obs
    time_obs = fT1_obs.variables["TIME"][:]*24 
    print ('time obs:',time_obs)
    ssh_obs = fT1_obs.variables["sossheig"][:,0]
    print(ssh_obs)
    # times
    te_obs = ssh_obs.shape[0]
    print ('te_obs=',te_obs)
    te_mod = ssh_mod.shape[0]
    print ('te_mod=',te_mod)
    ts_obs=0
    ts_mod=0
    # Fit
    print("Fit obs...")
    fitted_obs, cov_obs = curve_fit(octuple_obs,time_obs[ts_obs:480],ssh_obs[ts_obs:480])
    print ("Fit mod...")
    fitted_mod, cov_mod = curve_fit(octuple_mod,time_mod[ts_mod:480],ssh_mod[ts_mod:480])

    txtname=path+"pre_obs_mod"+stations_mod[stn]+".txt"
    np.savetxt(txtname, (fitted_obs,fitted_mod), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")
    #np.savetxt(txtname, (fitted_mod), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")   

    print ('OBS Val corr...')
    for idxA_fit in (0,2,4,6,8,10,12,14):
        if fitted_obs[idxA_fit] < 0:
          print (idxA_fit) 
          print ('<0')
          fitted_obs[idxA_fit]= -fitted_obs[idxA_fit]
          fitted_obs[idxA_fit+1]=fitted_obs[idxA_fit+1]+180 
    for idxP_fit in (1,3,5,7,9,11,13,15):    
         if fitted_obs[idxP_fit] < 0:
           print (idxP_fit)
           print ('<0') 
           fitted_obs[idxP_fit]=fitted_obs[idxP_fit]+360
         elif fitted_obs[idxP_fit] > 360:
           print (idxP_fit)
           print ('>360')
           fitted_obs[idxP_fit]=fitted_obs[idxP_fit]-360

    print ('MOD Val corr...')
    for idxA_fit in (0,2,4,6,8,10,12,14):
        if fitted_mod[idxA_fit] < 0:
           print (idxA_fit)
           print ('<0')
           fitted_mod[idxA_fit]= -fitted_mod[idxA_fit]
           fitted_mod[idxA_fit+1]=fitted_mod[idxA_fit+1]+180
    for idxP_fit in (1,3,5,7,9,11,13,15):
        if fitted_mod[idxP_fit] < 0:
           print (idxP_fit)
           print ('<0')
           fitted_mod[idxP_fit]=fitted_mod[idxP_fit]+360
        elif fitted_mod[idxP_fit] > 360:
           print (idxP_fit)
           print ('>360')
           fitted_mod[idxP_fit]=fitted_mod[idxP_fit]-360


    txtname=path+"obs_mod"+stations_mod[stn]+".txt"
    np.savetxt(txtname, (fitted_obs,fitted_mod), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")
    #np.savetxt(txtname, (fitted_mod), fmt='%5.3f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")

       
######    if fitted_obs[0] < 0:
######        fitted_obs[0] = -fitted_obs[0]
######        #fitted_obs[1] = fitted_obs[1]+180
######    
######    # Add nodal corrections:
######    M2_amp_obs.append(fitted_obs[0]*M2ft)
######    pha_obs = fitted_obs[1]+M2uvt
######    # Phase in [0,360]deg:
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    elif pha_obs < 0:
######        pha_obs = pha_obs+360
######    M2_pha_obs.append(pha_obs)
######    
######    if fitted_obs[2] < 0:
######        fitted_obs[2] = - fitted_obs[2]
######        fitted_obs[3] = fitted_obs[3] + 180
######    K1_amp_obs.append(fitted_obs[2]*K1ft)
######    pha_obs = fitted_obs[3] + K1uvt
######    if  pha_obs > 360:
######        pha_obs = pha_obs-360
######    K1_pha_obs.append(pha_obs)  
######    
######    if fitted_obs[4] < 0:
######        fitted_obs[4] = -fitted_obs[4]
######        fitted_obs[5] = fitted_obs[5]+180
######    O1_amp_obs.append(fitted_obs[4]*O1ft)
######    pha_obs= fitted_obs[5]+O1uvt
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    O1_pha_obs.append(pha_obs) 
######    
######    if fitted_obs[6] < 0:
######        fitted_obs[6] = -fitted_obs[6]
######        fitted_obs[7] = fitted_obs[7]+180
######    S2_amp_obs.append(fitted_obs[6]*S2ft)
######    pha_obs= fitted_obs[7]+S2uvt
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    S2_pha_obs.append(pha_obs) 
######    
######    
######    #####
######    if fitted_obs[8] < 0:
######            fitted_obs[8] = -fitted_obs[8]
######            fitted_obs[9] = fitted_obs[9]+180
######    P1_amp_obs.append(fitted_obs[8]*P1ft)
######    pha_obs= fitted_obs[9]+P1uvt
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    P1_pha_obs.append(pha_obs) 
######    
######    if fitted_obs[10] < 0:
######            fitted_obs[10] = -fitted_obs[10]
######            fitted_obs[11] = fitted_obs[11]+180
######    N2_amp_obs.append(fitted_obs[10]*N2ft)
######    pha_obs= fitted_obs[11]+N2uvt
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    N2_pha_obs.append(pha_obs) 
######    
######    if fitted_obs[12] < 0:
######            fitted_obs[12] = -fitted_obs[12]
######            fitted_obs[13] = fitted_obs[13]+180
######    Q1_amp_obs.append(fitted_obs[12]*Q1ft)
######    pha_obs= fitted_obs[13]+Q1uvt
######    if  pha_obs > 360:
######        pha_obs=pha_obs-360
######    Q1_pha_obs.append(pha_obs) 
######    
######    if fitted_obs[14] < 0:
######            fitted_obs[14] = -fitted_obs[14]
######            fitted_obs[15] = fitted_obs[15]+180
######    K2_amp_obs.append(fitted_obs[14]*K2ft)
######    pha_obs= fitted_obs[15]+K2uvt
######    if pha_obs > 360:
######        pha_obs = pha_obs-360
######    if pha_obs < 0:
######        pha_obs = pha_obs + 360
######    K2_pha_obs.append(pha_obs) 
######    
######    print('##########FIT OBS')
######    print(fitted_obs[:])
######    #file.write('# FIT OBS')
######    #np.savetxt('file_obsmod.txt', zip(fitted_obs,fitted_mod), fmt="%5.2f %5.2f")
######    
######    if fitted_mod[0] < 0:
######        fitted_mod[0] = -fitted_mod[0]
######        #fitted_mod[1] = fitted_mod[1]+180
######    
######    M2_amp_mod.append(fitted_mod[0]*M2ft)
######    pha_mod = fitted_mod[1]+M2uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    elif pha_mod < 0:
######        pha_mod = pha_mod+360
######    M2_pha_mod.append(pha_mod)
######    
######    if fitted_mod[2] < 0:
######        fitted_mod[2] = - fitted_mod[2]
######        fitted_mod[3] = fitted_mod[3] + 180
######    K1_amp_mod.append(fitted_mod[2]*K1ft)
######    pha_mod = fitted_mod[3] + K1uvt
######    if  pha_mod > 360:
######        pha_mod = pha_mod-360
######    K1_pha_mod.append(pha_mod)
######    
######    if fitted_mod[4] < 0:
######        fitted_mod[4] = -fitted_mod[4]
######        fitted_mod[5] = fitted_mod[5]+180
######    O1_amp_mod.append(fitted_mod[4]*O1ft)
######    pha_mod= fitted_mod[5]+O1uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    O1_pha_mod.append(pha_mod)
######    
######    if fitted_mod[6] < 0:
######        fitted_mod[6] = -fitted_mod[6]
######        fitted_mod[7] = fitted_mod[7]+180
######    S2_amp_mod.append(fitted_mod[6]*S2ft)
######    pha_mod= fitted_mod[7]+S2uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    S2_pha_mod.append(pha_mod)
######    
######    
######    #####
######    if fitted_mod[8] < 0:
######            fitted_mod[8] = -fitted_mod[8]
######            fitted_mod[9] = fitted_mod[9]+180
######    P1_amp_mod.append(fitted_mod[8]*P1ft)
######    pha_mod= fitted_mod[9]+P1uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    P1_pha_mod.append(pha_mod)
######    
######    if fitted_mod[10] < 0:
######            fitted_mod[10] = -fitted_mod[10]
######            fitted_mod[11] = fitted_mod[11]+180
######    N2_amp_mod.append(fitted_mod[10]*N2ft)
######    pha_mod= fitted_mod[11]+N2uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    N2_pha_mod.append(pha_mod)
######    
######    if fitted_mod[12] < 0:
######            fitted_mod[12] = -fitted_mod[12]
######            fitted_mod[13] = fitted_mod[13]+180
######    Q1_amp_mod.append(fitted_mod[12]*Q1ft)
######    pha_mod= fitted_mod[13]+Q1uvt
######    if  pha_mod > 360:
######        pha_mod=pha_mod-360
######    Q1_pha_mod.append(pha_mod)
######    
######    if fitted_mod[14] < 0:
######            fitted_mod[14] = -fitted_mod[14]
######            fitted_mod[15] = fitted_mod[15]+180
######    K2_amp_mod.append(fitted_mod[14]*K2ft)
######    pha_mod= fitted_mod[15]+K2uvt
######    if pha_mod > 360:
######        pha_mod = pha_mod-360
######    if pha_mod < 0:
######        pha_mod = pha_mod + 360
######    K2_pha_mod.append(pha_mod)
######    
######    print('##########FIT MOD')
######    print(fitted_mod[:])
######
######    txtname=path+"obs_mod"+stations_mod[stn]+".txt"
######    # Write arrays to file
######    #obs_arr=M2_amp_obs+M2_pha_obs+K1_amp_obs+K1_pha_obs+O1_amp_obs+O1_pha_obs+S2_amp_obs+S2_pha_obs+P1_amp_obs+P1_pha_obs+N2_amp_obs+N2_pha_obs+Q1_amp_obs+Q1_pha_obs+K2_amp_obs+K2_pha_obs
######    #mod_arr=M2_amp_mod+M2_pha_mod+K1_amp_mod+K1_pha_mod+O1_amp_mod+O1_pha_mod+S2_amp_mod+S2_pha_mod+P1_amp_mod+P1_pha_mod+N2_amp_mod+N2_pha_mod+Q1_amp_mod+Q1_pha_mod+K2_amp_mod+K2_pha_mod
######    #print(obs_arr)
######    #obs_arr.shape
######    #np.squeeze(mod_arr).shape
######    ##print(obs_arr)
######    np.savetxt(txtname, (fitted_obs,fitted_mod), fmt='%5.2f', delimiter=" ", header=" M2amp M2pha K1amp K1pha O1amp O1pha S2amp S2pha P1amp P1pha N2amp N2pha Q1amp Q1pha K2amp K2pha", comments="#")
######
#######first_col=['M2amp''M2pha''K1amp''K1pha''O1amp''O1pha''S2amp''S2pha''P1amp''P1pha''N2amp''N2pha''Q1amp''Q1pha''K2amp''K2pha']
#######M2_amp_obs,M2_pha_obs,K1_amp_obs,K1_pha_obs,O1_amp_obs,O1_pha_obs,S2_amp_obs,S2_pha_obs,P1_amp_obs,P1_pha_obs,N2_amp_obs,N2_pha_obs,Q1_amp_obs,Q1_pha_obs,K2_amp_obs,K2_pha_obs
#######M2_amp_mod,M2_pha_mod,K1_amp_mod,K1_pha_mod,O1_amp_mod,O1_pha_mod,S2_amp_mod,S2_pha_mod,P1_amp_mod,P1_pha_mod,N2_amp_mod,N2_pha_mod,Q1_amp_mod,Q1_pha_mod,K2_amp_mod,K2_pha_mod
######
#########
####### PLOT
####### M2 lin reg
#########labels=['Barcelona']
#########split1=8; split2=20
#########fig=plt.scatter(M2_amp_mod,M2_amp_obs,'M2',figsize=(14,6),split1=split1,split2=split2, labels=labels)
######
######
######
######
######
######
######
#######,M2_pha_mod,M2_pha_obs,'M2',figsize=(14,6),
#######                                  split1=split1,split2=split2, labels=labels)
######
#######ax_amp,ax_pha = fig.axes
#######min_value, max_value = ax_amp.set_xlim(0, 1.2)
#######ax_amp.plot([min_value, max_value], [min_value, max_value], color='red',lw=2)
######
#######min_value, max_value = ax_pha.set_xlim(0, 360)
#######ax_pha.plot([min_value, max_value], [min_value, max_value], color='red',lw=2)
######
#######plt.savefig('M2_Pha_Amp.jpg'
######
######
######
######
######tidal_c=['M2','N2','K2','S2','K1','O1','Q1','P1']
######
#######for comp in (0,1,2,3,4,5,6,7):
#######    print(tidal_c[comp])
######
######plt.figure(figsize=(10,6))
#########
######plt.subplot(1,2,1)
######plt.plot(np.array(M2_amp_mod)*1.03, '-bo', label = 'model')
######plt.plot(np.array(M2_amp_obs)*1.03, 'r-o', label = 'observation')
######plt.title('M2 Amplitude [m]')
######plt.legend( loc='upper left' )
#########
#########
######plt.subplot(1,2,2)
######plt.plot(np.array(M2_pha_mod)+2.3, '-bo', label = 'model')
######plt.plot(np.array(M2_pha_obs)+2.3, 'r-o', label = 'observation')
######plt.title('M2 Phase [deg]')
######
######plt.savefig(path+'M2_obsmod.jpg')
###
###
###plt.subplot(3,2,3)
#### use the un-wrap function to plot the M2 phase more smoothly
###pha_uwm = 180./np.pi prova2.py prova3.py prova.py tidetools_SSM.py np.unwrap(np.array(M2_pha)*np.pi/180.)
###plt.plot(pha_uwm+2.3, '-bo', label = 'model')
####pha_uw = 180./np.pi prova2.py prova3.py prova.py tidetools_SSM.py np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
####plt.plot(pha_uw, 'r-o', label = 'observation')
###plt.title('M2 Phase')
###
###plt.subplot(3,2,3)
#### use the un-wrap function to plot the M2 phase more smoothly
###pha_uwm = 180./np.pi prova2.py prova3.py prova.py tidetools_SSM.py np.unwrap(np.array(S2_pha)*np.pi/180.)
###plt.plot(pha_uwm+2.3, '-bo', label = 'model')
####pha_uw = 180./np.pi prova2.py prova3.py prova.py tidetools_SSM.py np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
####plt.plot(pha_uw, 'r-o', label = 'observation')
###plt.title('S2 Phase')
