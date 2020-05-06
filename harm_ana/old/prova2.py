#
# Script for DETIDING
# by Salish Sea MEOPAR
# https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis/raw/tip/compare_tides/Analysis8Components.ipynb
#

# imports
#%matplotlib inline
import matplotlib.pylab as plt
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

# pathname for data - all of the tide runs are stored in this directory
#path = '/data/nsoontie/MEOPAR/SalishSea/results/tides/'
path = '/work/ag15419/tmp/diaharm_ana/'

#the run we want to analyze
#runname = 'corr15'
runname = 'simt_ctrl0'

#joining the two string together
name = path +runname +'/'

print name

stations =  ['Gibraltar']

numsta=len(stations)

#constants and fitting
# M2
M2freq = 28.984106 # degrees per hour
M2freq = M2freq*np.pi/180. # radians per hour
#K1
K1freq = 15.041069*np.pi/180.
#O1
O1freq = 13.943036*np.pi/180.
#S2
S2freq = 30.000002*np.pi/180.
#P1
P1freq = 14.958932*np.pi/180.
#N2
N2freq = 28.439730*np.pi/180.
#Q1
Q1freq = 13.398661*np.pi/180.
#K2
K2freq = 30.082138*np.pi/180.

# initial phase calculation
M2ft = 1.036581
M2uvt = 0.009395
S2ft = 1.000000  
S2uvt = 0.000000

def double(x, M2amp, M2pha, S2amp, S2pha):
    return (M2amp*np.cos(M2freq*x-M2pha*np.pi/180.)+
            K1amp*np.cos(S2freq*x-S2pha*np.pi/180.))

fig, ax = plt.subplots(1,1,figsize=(12,5))
for stn in (0,4,14,23):
    print stations[stn]
    fT1 = NC.Dataset(name+stations[stn]+'.nc','r')
    time = fT1.variables["time_counter"][:]/3600.  # want hours not seconds
    ssh = fT1.variables["sossheig"][:,0,0]
    ax.plot(time,ssh)

print ssh.shape

M2_amp=[]; M2_pha=[]; S2_amp=[]; S2_pha=[]

ts = 150
te = ssh.shape[0]

for stn in range(numsta):
    fT1 = NC.Dataset(name+stations[stn]+'.nc','r')
    time = fT1.variables["time_counter"][:]/3600.  # want hours not seconds
    ssh = fT1.variables["sossheig"][:,0,0]

    ##fitted, cov = curve_fit(octuple,time[ts:te],ssh[ts:te]) 
    fitted, cov = curve_fit(double,time[ts:te],ssh[ts:te]) 
    if fitted[0] < 0:
        fitted[0] = -fitted[0]
        fitted[1] = fitted[1]+180

    M2_amp.append(fitted[0]*M2ft)
    pha = fitted[1]+M2uvt
    if  pha > 360:
        pha=pha-360
    elif pha < 0:
        pha = pha+360
    if stn == 6:
        print pha
    M2_pha.append(pha)

    if fitted[2] < 0:
        fitted[2] = -fitted[2]
        fitted[3] = fitted[3]+180
    S2_amp.append(fitted[2]*S2ft)
    pha= fitted[3]+S2uvt
    if  pha > 360:
        pha=pha-360
    S2_pha.append(pha) 

##
##def complex_diff(Ao,go,Am,gm):
##    #calculates complex differences between observations and model
##    #Ao, go - amplitude and phase from observations
##    #Am, gm - amplitude and phase from model
##    D = np.sqrt((Ao*np.cos(np.pi*go/180)-Am*np.cos(np.pi*gm/180))**2 + 
##                (Ao*np.sin(np.pi*go/180)-Am*np.sin(np.pi*gm/180))**2)
##    
##    return D
##
###R
##R_M2 = M2_amp/M2_amp_obs
##R_K1 = K1_amp/K1_amp_obs
###delta phi (adjust so between -180, 180)
##Dphi_M2 = M2_pha-M2_pha_obs; 
##Dphi_M2 = Dphi_M2 -360*(Dphi_M2>180) + 360*(Dphi_M2<-180)
##Dphi_K1 = K1_pha-K1_pha_obs
##Dphi_K1 = Dphi_K1 -360*(Dphi_K1>180) + 360*(Dphi_K1<-180)
###Complex differences
##D_M2= complex_diff(np.array(M2_amp_obs),np.array(M2_pha_obs), np.array(M2_amp)*1.03,np.array(M2_pha)+2.3)
##D_K1= complex_diff(np.array(K1_amp_obs),np.array(K1_pha_obs), np.array(K1_amp)*0.99,np.array(K1_pha)-0.5)
##D_O1= complex_diff(np.ma.array(O1_amp_obs),np.ma.array(O1_pha_obs), np.ma.array(O1_amp),np.ma.array(O1_pha))
##D_S2= complex_diff(np.ma.array(S2_amp_obs),np.ma.array(S2_pha_obs), np.ma.array(S2_amp),np.ma.array(S2_pha))
##D_P1= complex_diff(np.ma.array(P1_amp_obs),np.ma.array(P1_pha_obs), np.ma.array(P1_amp),np.ma.array(P1_pha))
##D_N2= complex_diff(np.ma.array(N2_amp_obs),np.ma.array(N2_pha_obs), np.ma.array(N2_amp),np.ma.array(N2_pha))
##D_Q1= complex_diff(np.ma.array(Q1_amp_obs),np.ma.array(Q1_pha_obs), np.ma.array(Q1_amp),np.ma.array(Q1_pha))
##D_K2= complex_diff(np.ma.array(K2_amp_obs),np.ma.array(K2_pha_obs), np.ma.array(K2_amp),np.ma.array(K2_pha))
##
##outfile = runname+'.csv'
##
##with open(outfile, 'wb') as csvfile:
##    writer = csv.writer(csvfile, delimiter=',')
##    writer.writerow([
##            'Station Name', 
##            'R (M2)', 'Delta phi (M2)', 'D (M2)',
##            'R (K1)', 'Delta phi (K1)', 'D (K1)'
##        ])
##    for stn in range(numsta):
##        location = stations_obs[stn]
##        writer.writerow([stations_obs[stn],
##                        R_M2[stn], Dphi_M2[stn], D_M2[stn],
##                        R_K1[stn], Dphi_K1[stn], D_K1[stn]])
##
##    #write averages and rms
##    writer.writerow(['Mean Difference',
##                    mean(M2_amp-M2_amp_obs),mean(Dphi_M2),mean(D_M2), 
##                    mean(K1_amp-K1_amp_obs),mean(Dphi_K1),mean(D_K1)])
##    writer.writerow(['RMS Difference',
##                    rms(M2_amp-M2_amp_obs),rms(Dphi_M2),rms(D_M2), 
##                    rms(K1_amp-K1_amp_obs),rms(Dphi_K1),rms(D_K1)])
##    #without the north
##    writer.writerow(['Mean Difference no North no PR',
##                    mean(M2_amp[1:split2]-M2_amp_obs[1:split2]),mean(Dphi_M2[1:split2]),mean(D_M2[1:split2]), 
##                    mean(K1_amp[1:split2]-K1_amp_obs[1:split2]),mean(Dphi_K1[1:split2]),mean(D_K1[1:split2])])
##    writer.writerow(['RMS Difference no North no PR',
##                    rms(M2_amp[1:split2]-M2_amp_obs[1:split2]),rms(Dphi_M2[1:split2]),rms(D_M2[1:split2]), 
##                    rms(K1_amp[1:split2]-K1_amp_obs[1:split2]),rms(Dphi_K1[1:split2]),rms(D_K1[1:split2])])

plt.figure(figsize=(20,12))

plt.subplot(3,2,1)
plt.plot(np.array(M2_amp)*1.03, '-bo', label = 'model')
#plt.plot(M2_amp_obs, 'r-o', label = 'observation')
plt.title('M2 Amplitude')
plt.legend( loc='upper left' )

##plt.subplot(3,2,2)
##plt.plot(np.array(K1_amp)*0.99, '-bo', label = 'model')
##plt.plot(K1_amp_obs, 'r-o', label = 'observation')
##plt.title('K1 Amplitude')

plt.subplot(3,2,2)
plt.plot(np.array(S2_amp)*0.99, '-bo', label = 'model')
##plt.plot(K1_amp_obs, 'r-o', label = 'observation')
plt.title('S2 Amplitude')


plt.subplot(3,2,3)
# use the un-wrap function to plot the M2 phase more smoothly
pha_uwm = 180./np.pi * np.unwrap(np.array(M2_pha)*np.pi/180.)
plt.plot(pha_uwm+2.3, '-bo', label = 'model')
#pha_uw = 180./np.pi * np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
#plt.plot(pha_uw, 'r-o', label = 'observation')
plt.title('M2 Phase')

plt.subplot(3,2,3)
# use the un-wrap function to plot the M2 phase more smoothly
pha_uwm = 180./np.pi * np.unwrap(np.array(S2_pha)*np.pi/180.)
plt.plot(pha_uwm+2.3, '-bo', label = 'model')
#pha_uw = 180./np.pi * np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
#plt.plot(pha_uw, 'r-o', label = 'observation')
plt.title('S2 Phase')



##plt.subplot(3,2,4)
##pha_uw = 180./np.pi * np.unwrap(np.array(K1_pha)*np.pi/180.)
##plt.plot(pha_uw-0.5, '-bo', label = 'model')
##plt.plot(K1_pha_obs, 'r-o', label = 'observation')
##plt.title('K1 Phase')

##plt.subplot(3,2,5)
##plt.plot(D_M2, '-bo', label = 'M2')
##plt.plot(D_K1, '-go', label = 'K1')
##plt.plot((0,30),(0.05,0.05),'k')
##plt.plot((0,30),(0.10,0.10),'r')
##plt.title('D error')
##plt.legend( loc='upper left' )


##M2_amp_topogbottfric = M2_amp
##M2_pha_topogbottfric = pha_uwm
##
##M2_amp_rnshlat2 = M2_amp
##M2_pha_rnshlat2 = pha_uwm
##
##M2_amp_bot1em2B = M2_amp
##M2_pha_bot1em2B = pha_uwm
##
##M2_amp_bot1em2 = M2_amp
##M2_pha_bot1em2 = pha_uwm
##
##M2_amp_corr15 = M2_amp
##M2_pha_corr15 = pha_uwm
##
##M2_amp_topog = M2_amp
##M2_pha_topog = pha_uwm
##
##M2_amp_rnshlat = M2_amp
##M2_pha_rnshlat = pha_uwm
##
##plt.figure(figsize=(12,5))
##plt.plot(np.array(M2_amp_topog)*0.97, '-bo', label = 'topog')
##plt.plot(M2_amp_obs, 'r-s', label = 'observation')
##plt.plot(np.array(M2_amp_rnshlat2)*1.08, '-m^')
##plt.plot(np.array(M2_amp_bot1em2B)*1.07, '-c^')
##plt.plot(M2_amp_corr15, '-go', label='corr15')
##plt.plot(np.array(M2_amp_topogbottfric)*1.03, '-yo')
##makeit = np.array(M2_amp_corr15) + (np.array(M2_amp_bot1em2B)*1.07-np.array(M2_amp_corr15)) + (np.array(M2_amp_topog)*0.97-np.array(M2_amp_corr15))
##plt.plot(makeit, 'k*-')
##
##plt.figure(figsize=(12,5))
##plt.plot(np.array(M2_pha_topog)+2.9, '-bo', label = 'topog')
##pha_uw = 180./np.pi * np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
##plt.plot(pha_uw, 'r-s', label = 'observation')
##plt.plot(np.array(M2_pha_rnshlat2)-0.8, '-m^')
##plt.plot(np.array(M2_pha_bot1em2B)-0.5, '-c^')
##plt.plot(M2_pha_corr15, '-go', label = 'corr15')
##plt.plot(np.array(M2_pha_topogbottfric)+2, '-yo')
##makeit = np.array(M2_pha_corr15) + (np.array(M2_pha_bot1em2B)-0.8-np.array(M2_pha_corr15)) + (np.array(M2_pha_topog)+2.9-np.array(M2_pha_corr15))
##plt.plot(makeit, 'k*-')
##
##cmap = plt.get_cmap('PuBu')
##cmap.set_bad('burlywood')
##
##fig,axs=plt.subplots(4, 2, figsize=(8,20))
##
##constituent = ('M2', 'K1', 'O1', 'S2', 'P1', 'N2', 'Q1', 'K2')
##error_D = (D_M2, D_K1, D_O1, D_S2, D_P1, D_N2, D_Q1, D_K2)
##
##
##for row in range(4):
##
##    for ax, error_D1, const in zip(axs[row], error_D[row*2:row*2+2], constituent[row*2:row*2+2]):
##        ax.pcolormesh(X,Y,bathy,cmap='PuBu')
##
##        for stn in range(numsta):
##            location = stations_obs[stn]
##            lon=-harm_obs.lon[harm_obs.site==location]
##            lat=harm_obs.lat[harm_obs.site==location]
##            if error_D1 [stn] <= 0.05:
##                ax.plot(lon,lat,'og',label=location,markersize=10,markeredgecolor='g')
##            if error_D1 [stn] > 0.1:
##                ax.plot(lon,lat,'or',label=location,markersize=10,markeredgecolor='r')
##            if 0.1 >= error_D1[stn] > 0.05:
##                ax.plot(lon,lat,'oy',label=location,markersize=10,markeredgecolor='y')
##        
##            ax.annotate(stn, xy = (lon,lat), xytext = (5,5),ha = 'right', va = 'bottom',
##                textcoords = 'offset points')
##            ax.set_title(const)
##        ax.axis([-126.1,-122,47,51])
##   
##
##
##
##fig, axs = plt.subplots(6,2,figsize=(10,15))
##axs[0,0].plot(np.array(O1_amp)/np.array(K1_amp), '-bo', label = 'model')
##axs[0,0].plot((0,28),(0.560,0.560), 'r-', label = 'observation')
##axs[0,0].set_title('O1/K1 Amplitude')
##axs[0,1].plot(np.array(O1_pha)-np.array(K1_pha), '-bo', label = 'model')
##axs[0,1].plot((0,28),(-22.9,-22.9), 'r-', label = 'observation')
##axs[0,1].set_title('O1-K1 Phase')
##
##axs[1,0].plot(np.array(S2_amp)/np.array(M2_amp), '-bo', label = 'model')
##axs[1,0].plot((0,28),(0.249,0.249), 'r-', label = 'observation')
##axs[1,0].set_title('S2/M2 Amplitude')
##pha_uw = 180./np.pi * np.unwrap((np.array(S2_pha)-np.array(M2_pha))*np.pi/180.)
##axs[1,1].plot(pha_uw, '-bo', label = 'model')
##axs[1,1].plot((0,28),( 28.7, 28.7), 'r-', label = 'observation')
##axs[1,1].set_title('S2-M2 Phase')
##
##axs[2,0].plot(np.array(P1_amp)/np.array(K1_amp), '-bo', label = 'model')
##axs[2,0].plot((0,28),(0.311,0.311), 'r-', label = 'observation')
##axs[2,0].set_title('P1/K1 Amplitude')
##
##pha_uw = 180./np.pi * np.unwrap((np.array(P1_pha)-np.array(K1_pha))*np.pi/180.)
##axs[2,1].plot(pha_uw, '-bo', label = 'model')
##axs[2,1].plot((0,28),(-3,-3), 'r-', label = 'observation')
##axs[2,1].set_title('P1-K1 Phase')
##
##axs[3,0].plot(np.array(N2_amp)/np.array(M2_amp), '-bo', label = 'model')
##axs[3,0].plot((0,28),(0.200,0.200), 'r-', label = 'observation')
##axs[3,0].set_title('N2/M2 Amplitude')
##
##pha_uw = 180./np.pi * np.unwrap((np.array(N2_pha)-np.array(M2_pha))*np.pi/180.)
##axs[3,1].plot(pha_uw, '-bo', label = 'model')
##axs[3,1].plot((0,28),(-28.3, -28.3), 'r-', label = 'observation')
##axs[3,1].set_title('N2-M2 Phase')
##
##axs[4,0].plot(np.array(Q1_amp)/np.array(K1_amp), '-bo', label = 'model')
##axs[4,0].plot((0,28),(0.089,0.089), 'r-', label = 'observation')
##axs[4,0].set_title('Q1/K1 Amplitude')
##
##pha_uw = 180./np.pi * np.unwrap((np.array(Q1_pha)-np.array(K1_pha))*np.pi/180.)
##axs[4,1].plot(pha_uw+360., '-bo', label = 'model')
##axs[4,1].plot((0,28),(-27.3,-27.3), 'r-', label = 'observation')
##axs[4,1].set_title('Q1-K1 Phase')
##
##axs[5,0].plot(np.array(K2_amp)/np.array(M2_amp), '-bo', label = 'model')
##axs[5,0].plot((0,28),(0.068,0.068), 'r-', label = 'observation')
##axs[5,0].set_title('K2/M2 Amplitude')
##
##pha_uw = 180./np.pi * np.unwrap((np.array(K2_pha)-np.array(M2_pha))*np.pi/180.)
##axs[5,1].plot(pha_uw, '-bo', label = 'model')
##axs[5,1].plot((0,28),(28.7, 28.7), 'r-', label = 'observation')
##axs[5,1].set_title('K2-M2 Phase')
##
##print te
##
##sample = 17
##start = np.zeros(sample)
##tend = np.zeros(sample)
##for i in range(sample):
##    start[i] = 196+(480-196)*np.random.rand()
##    tend[i] = te-(480-196)*np.random.rand()
##print start
##print tend
##timelength = (tend-start)/96.
##print np.mean(timelength),2*np.std(timelength)
##print time[start[1]:tend[1]]
##
##
###allocate space for our arrays
##M2_amp=np.zeros((numsta,sample)); M2_pha=np.zeros((numsta,sample))
##K1_amp=np.zeros((numsta,sample)); K1_pha=np.zeros((numsta,sample))
##O1_amp=np.zeros((numsta,sample)); O1_pha=np.zeros((numsta,sample))
##S2_amp=np.zeros((numsta,sample)); S2_pha=np.zeros((numsta,sample))
##P1_amp=np.zeros((numsta,sample)); P1_pha=np.zeros((numsta,sample))
##N2_amp=np.zeros((numsta,sample)); N2_pha=np.zeros((numsta,sample))
##Q1_amp=np.zeros((numsta,sample)); Q1_pha=np.zeros((numsta,sample))
##K2_amp=np.zeros((numsta,sample)); K2_pha=np.zeros((numsta,sample))
##
##
##
##for it,tst,tet in zip(range(sample),start,tend):
##  
##    for stn in range(numsta):
##        fT1 = NC.Dataset(name+stations[stn]+'.nc','r')
##        time = fT1.variables["time_counter"][:]/3600.  # want hours not seconds
##        ssh = fT1.variables["sossheig"][:,0,0]
##
##        fitted, cov = curve_fit(octuple,time[tst:tet],ssh[tst:tet]) 
##        if fitted[0] < 0:
##            fitted[0] = -fitted[0]
##            fitted[1] = fitted[1]+180
##
##        M2_amp[stn,it] = fitted[0]*M2ft
##        pha = fitted[1]+M2uvt
##        if  pha > 360:
##            pha=pha-360
##        M2_pha[stn,it] = pha
##        
##        if fitted[2] < 0:
##            fitted[2] = -fitted[2]
##            fitted[3] = fitted[3]+180
##
##        K1_amp[stn,it] = fitted[2]*K1ft
##        pha= fitted[3]+K1uvt
##        if  pha > 360:
##            pha=pha-360
##        K1_pha[stn,it]= pha   
##                
##        if fitted[4] < 0:
##            fitted[4] = -fitted[4]
##            fitted[5] = fitted[5]+180
##        O1_amp[stn,it] =fitted[4]*O1ft
##        pha= fitted[5]+O1uvt
##        if  pha > 360:
##            pha=pha-360
##        O1_pha[stn,it]= pha 
##        
##        if fitted[6] < 0:
##            fitted[6] = -fitted[6]
##            fitted[7] = fitted[7]+180
##        S2_amp[stn,it] =fitted[6]*S2ft
##        pha= fitted[7]+S2uvt
##        if  pha > 360:
##            pha=pha-360
##        S2_pha[stn,it]= pha 
##        
##        if fitted[8] < 0:
##            fitted[8] = -fitted[8]
##            fitted[9] = fitted[9]+180
##        P1_amp[stn,it] = fitted[8]*P1ft
##        pha= fitted[9]+P1uvt
##        if  pha > 360:
##            pha=pha-360
##        P1_pha[stn,it] =pha 
##    
##        if fitted[10] < 0:
##            fitted[10] = -fitted[10]
##            fitted[11] = fitted[11]+180
##        N2_amp[stn,it] = fitted[10]*N2ft
##        pha= fitted[11]+N2uvt
##        if  pha > 360:
##            pha=pha-360
##        N2_pha[stn,it] = pha
##        
##        if fitted[12] < 0:
##            fitted[12] = -fitted[12]
##            fitted[13] = fitted[13]+180
##        Q1_amp[stn,it] = fitted[12]*Q1ft
##        pha= fitted[13]+Q1uvt
##        if  pha > 360:
##            pha=pha-360
##        Q1_pha[stn,it] = pha
##        
##        if fitted[14] < 0:
##            fitted[14] = -fitted[14]
##            fitted[15] = fitted[15]+180
##        K2_amp[stn,it] = fitted[14]*K2ft
##        pha= fitted[15]+K2uvt
##        if  pha > 360:
##            pha=pha-360
##        K2_pha[stn,it] = pha
##    
##
##jdef = range(3)
##south = range(14,18)
##north = range(29,31)
##print 'M2'
##print '     JdeFuca'
##print np.mean(M2_amp[jdef]),2*np.std(np.mean(M2_amp[jdef],axis=0))
##print np.mean(M2_amp_obs[jdef]), np.mean(M2_amp_obs[jdef])-np.mean(M2_amp[jdef])
##print np.mean(M2_amp_obs[jdef])/np.mean(M2_amp[jdef])
##print np.mean(M2_pha[jdef]),2*np.std(np.mean(M2_pha[jdef],axis=0))
##print np.mean(M2_pha_obs[jdef]), np.mean(M2_pha_obs[jdef])-np.mean(M2_pha[jdef])
##print '     South'
##print np.mean(M2_amp[south]),2*np.std(np.mean(M2_amp[south],axis=0))
##print np.mean(M2_amp_obs[south]), np.mean(M2_amp_obs[south])-np.mean(M2_amp[south])
##print np.mean(M2_amp_obs[south])/np.mean(M2_amp[south])
##print np.mean(M2_pha[south]),2*np.std(np.mean(M2_pha[south],axis=0))
##print np.mean(M2_pha_obs[south]), np.mean(M2_pha_obs[south])-np.mean(M2_pha[south])
##print '     North'
##print np.mean(M2_amp[north]),2*np.std(np.mean(M2_amp[north],axis=0))
##print np.mean(M2_amp_obs[north]), np.mean(M2_amp_obs[north])-np.mean(M2_amp[north])
##print np.mean(M2_amp_obs[north])/np.mean(M2_amp[north])
##print np.mean(M2_pha[north]),2*np.std(np.mean(M2_pha[north],axis=0))
##print np.mean(M2_pha_obs[north]), np.mean(M2_pha_obs[north])-np.mean(M2_pha[north])
##print '==============================================='
##print 'K1'
##print '     JdeFuca'
##print np.mean(K1_amp[jdef]),2*np.std(np.mean(K1_amp[jdef],axis=0))
##print np.mean(K1_amp_obs[jdef]), np.mean(K1_amp_obs[jdef])-np.mean(K1_amp[jdef])
##print np.mean(K1_amp_obs[jdef])/np.mean(K1_amp[jdef])
##print np.mean(K1_pha[jdef]),2*np.std(np.mean(K1_pha[jdef],axis=0))
##print np.mean(K1_pha_obs[jdef]), np.mean(K1_pha_obs[jdef])-np.mean(K1_pha[jdef])
##print '     South'
##print np.mean(K1_amp[south]),2*np.std(np.mean(K1_amp[south],axis=0))
##print np.mean(K1_amp_obs[south]), np.mean(K1_amp_obs[south])-np.mean(K1_amp[south])
##print np.mean(K1_amp_obs[south])/np.mean(K1_amp[south])
##print np.mean(K1_pha[south]),2*np.std(np.mean(K1_pha[south],axis=0))
##print np.mean(K1_pha_obs[south]), np.mean(K1_pha_obs[south])-np.mean(K1_pha[south])
##print '     North'
##print np.mean(K1_amp[north]),2*np.std(np.mean(K1_amp[north],axis=0))
##print np.mean(K1_amp_obs[north]), np.mean(K1_amp_obs[north])-np.mean(K1_amp[north])
##print np.mean(K1_amp_obs[north])/np.mean(K1_amp[north])
##print np.mean(K1_pha[north]),2*np.std(np.mean(K1_pha[north],axis=0))
##print np.mean(K1_pha_obs[north]), np.mean(K1_pha_obs[north])-np.mean(K1_pha[north])
##print '==============================================='
##
##print 'O1'
##print '     South'
##print np.mean(O1_amp[south]/K1_amp[south]),2*np.std(np.mean(O1_amp[south]/K1_amp[south],axis=0))
##print np.mean(O1_amp_obs[south]/K1_amp_obs[south]), (np.mean(O1_amp_obs[south]/K1_amp_obs[south])
##                                                    -np.mean(O1_amp[south]/K1_amp[south]))
##print np.mean(O1_amp_obs[south]/K1_amp_obs[south])/np.mean(O1_amp[south]/K1_amp[south])
##print np.mean(O1_pha[south]-K1_pha[south]),2*np.std(np.mean(O1_pha[south]-K1_pha[south],axis=0))
##print np.mean(O1_pha_obs[south]-K1_pha_obs[south]), (np.mean(O1_pha_obs[south]-K1_pha_obs[south])
##                                                     -np.mean(O1_pha[south]-K1_pha[south]))
##print '     North'
##print np.mean(O1_amp[north]/K1_amp[north]),2*np.std(np.mean(O1_amp[north]/K1_amp[north],axis=0))
##print np.mean(O1_amp_obs[north]/K1_amp_obs[north]), (np.mean(O1_amp_obs[north]/K1_amp_obs[north])
##                                                    -np.mean(O1_amp[north]/K1_amp[north]))
##print np.mean(O1_amp_obs[north]/K1_amp_obs[north])/np.mean(O1_amp[north]/K1_amp[north])
##print np.mean(O1_pha[north]-K1_pha[north]),2*np.std(np.mean(O1_pha[north]-K1_pha[north],axis=0))
##print np.mean(O1_pha_obs[north]-K1_pha_obs[north]), (np.mean(O1_pha_obs[north]-K1_pha_obs[north])
##                                                     -np.mean(O1_pha[north]-K1_pha[north]))
##print '==============================================='
##print 'S2'
##code = ('south','north')
##for dir,dire in zip(code,(south,north)):
##    print dir
##    print np.mean(S2_amp[dire]/M2_amp[dire]),2*np.std(np.mean(S2_amp[dire]/M2_amp[dire],axis=0))
##    print np.mean(S2_amp_obs[dire]/M2_amp_obs[dire]), (np.mean(S2_amp_obs[dire]/M2_amp_obs[dire])
##                                                    -np.mean(S2_amp[dire]/M2_amp[dire]))
##    print np.mean(S2_amp_obs[dire]/M2_amp_obs[dire])/np.mean(S2_amp[dire]/M2_amp[dire])
##    unwrap = np.unwrap(np.array(S2_pha)*np.pi/180.)*180./np.pi
##    M2_un = np.unwrap(np.array(M2_pha)*np.pi/180.)*180./np.pi
##    plt.plot (unwrap[dire],'r',M2_un[dire],'b')
##    print np.mean(unwrap[dire]-M2_un[dire])+360.,2*np.std(np.mean(unwrap[dire]-M2_un[dire],axis=0))
##    print np.mean(S2_pha_obs[dire]-M2_pha_obs[dire]), (np.mean(S2_pha_obs[dire]-M2_pha_obs[dire])
##                                                     -np.mean(unwrap[dire]-M2_un[dire]))-360.
##
##const = ('P1', 'Q1')
##model_amp = (P1_amp, Q1_amp)
##model_pha = ()
##for const, model_amp, model_pha, obs_amp, obs_pha in zip(('P1','Q1'),
##                                                    (P1_amp, Q1_amp),(P1_pha, Q1_pha), 
##                                                  (P1_amp_obs, Q1_amp_obs), (P1_pha_obs, Q1_pha_obs)):
##    print const
##    for dir,dire in zip(code,(south,north)):
##        print dir
##        print np.mean(model_amp[dire]/K1_amp[dire]),2*np.std(np.mean(model_amp[dire]/K1_amp[dire],axis=0))
##        print np.mean(obs_amp[dire]/K1_amp_obs[dire]), (np.mean(obs_amp[dire]/K1_amp_obs[dire])
##                                                    -np.mean(model_amp[dire]/K1_amp[dire]))
##        print np.mean(obs_amp[dire]/K1_amp_obs[dire])/np.mean(model_amp[dire]/K1_amp[dire])
##        unwrap = np.unwrap(np.array(model_pha)*np.pi/180.)*180./np.pi
##        K1_un = np.unwrap(np.array(K1_pha)*np.pi/180.)*180./np.pi
##        print np.mean(unwrap[dire]-K1_un[dire]),2*np.std(np.mean(unwrap[dire]-K1_un[dire],axis=0))
##        print np.mean(obs_pha[dire]-K1_pha_obs[dire]), (np.mean(obs_pha[dire]-K1_pha_obs[dire])
##                                                     -np.mean(unwrap[dire]-K1_un[dire]))
##
##const = ('N2', 'K2')
##model_amp = (N2_amp, K2_amp)
##model_pha = ()
##for const, model_amp, model_pha, obs_amp, obs_pha in zip(('N2','K2'),
##                                                    (N2_amp, K2_amp),(N2_pha, K2_pha), 
##                                                  (N2_amp_obs, K2_amp_obs), (N2_pha_obs, K2_pha_obs)):
##    print const
##    for dir,dire in zip(code,(south,north)):
##        print dir
##        print np.mean(model_amp[dire]/M2_amp[dire]),2*np.std(np.mean(model_amp[dire]/M2_amp[dire],axis=0))
##        print np.mean(obs_amp[dire]/M2_amp_obs[dire]), (np.mean(obs_amp[dire]/M2_amp_obs[dire])
##                                                    -np.mean(model_amp[dire]/M2_amp[dire]))
##        print np.mean(obs_amp[dire]/M2_amp_obs[dire])/np.mean(model_amp[dire]/M2_amp[dire])
##        unwrap = model_pha#np.unwrap(np.array(model_pha)*np.pi/180.)*180./np.pi
##        M2_un = np.unwrap(np.array(M2_pha)*np.pi/180.)*180./np.pi
##        print unwrap[dire]
##        print np.mean(unwrap[dire]-M2_un[dire]),2*np.std(np.mean(unwrap[dire]-M2_un[dire],axis=0))
##        print np.mean(obs_pha[dire]-M2_pha_obs[dire]), (np.mean(obs_pha[dire]-M2_pha_obs[dire])
##                                                     -np.mean(unwrap[dire]-M2_un[dire]))
##
##
##
##fig,axs = plt.subplots(8,2,figsize=(15,25))
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(M2_pha[:,i])*np.pi/180.)
##    axs[0,1].plot(pha_uw ,'-ob', label = 'model')
##pha_uw = 180./np.pi * np.unwrap(np.array(M2_pha_obs)*np.pi/180.)
##axs[0,1].plot(pha_uw, 'r-*', label = 'observation')
##axs[0,1].set_title('M2 Phase')
##for i in range(sample):
##    axs[0,0].plot(M2_amp[:,i], '-bo', label = 'model')
##axs[0,0].plot(M2_amp_obs, 'r-*', label = 'observation')
##axs[0,0].set_title('M2 Amp')
##
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(K1_pha[:,i])*np.pi/180.)
##    axs[1,1].plot(pha_uw, '-bo', label = 'model')
##axs[1,1].plot(K1_pha_obs, 'r-*', label = 'observation')
##axs[1,1].set_title('K1 Phase')
##for i in range(sample):
##    axs[1,0].plot(K1_amp[:,i], '-bo', label = 'model')
##axs[1,0].plot(K1_amp_obs, 'r-*', label = 'observation')
##axs[1,0].set_title('K1 Amp')
##
##for i in range(sample):
##    if O1_pha[0,i] < 0:
##        O1_pha[0,i] = O1_pha[0,i] + 360
##    pha_uw = 180./np.pi * np.unwrap(np.array(O1_pha[:,i])*np.pi/180.)
##    axs[2,1].plot(pha_uw, '-bo', label = 'model')
##axs[2,1].plot(O1_pha_obs, 'r-*', label = 'observation', markersize = 15)
##axs[2,1].set_title('O1 Phase')
##for i in range(sample):
##    axs[2,0].plot(O1_amp[:,i], '-bo', label = 'model')
##axs[2,0].plot(O1_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[2,0].set_title('O1 Amp')
##
##for i in range(sample):
##    if S2_pha[0,i] < 0:
##        S2_pha[0,i] = S2_pha[0,i] + 360
##    pha_uw = 180./np.pi * np.unwrap(np.array(S2_pha[:,i])*np.pi/180.)
##    axs[3,1].plot(pha_uw, '-bo', label = 'model')
##pha_uw = 180./np.pi * np.unwrap(np.array(S2_pha_obs)*np.pi/180.)
##vsmall = 1e-6
##pha_uwm = np.ma.masked_array(pha_uw, mask=(abs(pha_uw-360)<vsmall))
##axs[3,1].plot(pha_uwm, 'r-*', label = 'observation', markersize = 15)
##axs[3,1].set_title('S2 Phase')
##for i in range(sample):
##    axs[3,0].plot(S2_amp[:,i], '-bo', label = 'model')
##axs[3,0].plot(S2_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[3,0].set_title('S2 Amp')
##
##for i in range(sample):
##    axs[4,0].plot(P1_amp[:,i], '-bo', label = 'model')
##axs[4,0].plot(P1_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[4,0].set_title('P1 Amp')
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(P1_pha[:,i])*np.pi/180.)
##    axs[4,1].plot(pha_uw, '-bo', label = 'model')
##axs[4,1].plot(P1_pha_obs, 'r-*', label = 'observation', markersize = 15)
##axs[4,1].set_title('P1 Phase')
##
##for i in range(sample):
##    axs[5,0].plot(N2_amp[:,i], '-bo', label = 'model')
##axs[5,0].plot(N2_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[5,0].set_title('N2 Amp')
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(N2_pha[:,i])*np.pi/180.)
##    axs[5,1].plot(pha_uw, '-bo', label = 'model')
##pha_uw = 180./np.pi * np.unwrap(np.array(N2_pha_obs)*np.pi/180.)
##pha_uwm = np.ma.masked_array(pha_uw, mask=(abs(pha_uw-360)<vsmall))
##axs[5,1].plot(pha_uwm, 'r-*', label = 'observation', markersize = 15)  
##axs[5,1].set_title('N2 Phase')
##
##for i in range(sample):
##    axs[6,0].plot(Q1_amp[:,i], '-bo', label = 'model')
##axs[6,0].plot(Q1_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[6,0].set_title('Q1 Amp')
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(Q1_pha[:,i])*np.pi/180.)
##    for j in range(numsta):
##        if pha_uw[j] < 0:
##            pha_uw[j] += 360
##    axs[6,1].plot(pha_uw, '-bo', label = 'model')
##axs[6,1].plot(Q1_pha_obs, 'r-*', label = 'observation', markersize = 15)
##axs[6,1].set_title('Q1 Phase')
##
##for i in range(sample):
##    axs[7,0].plot(K2_amp[:,i], '-bo', label = 'model')
##axs[7,0].plot(K2_amp_obs, 'r-*', label = 'observation', markersize = 15)
##axs[7,0].set_title('K2 Amp')
##for i in range(sample):
##    pha_uw = 180./np.pi * np.unwrap(np.array(K2_pha[:,i])*np.pi/180.)
##    for j in range(numsta):
##        if pha_uw[j] < 120:
##            pha_uw[j] += 360
##    axs[7,1].plot(pha_uw ,'-bo', label = 'model')
##for j in range(numsta):
##        if K2_pha_obs[j] < 120:
##            K2_pha_obs[j] += 360
##axs[7,1].plot(K2_pha_obs, 'r-*', label = 'observation', markersize = 15)
##axs[7,1].set_title('K2 Phase')
##
##fig,axs = plt.subplots(1,2,figsize=(15,10))
##K1mean=np.average(K1_pha, axis=1)
##K1max = np.max(K1_pha,axis=1)
##K1min = np.min(K1_pha,axis=1)
##asymmetric_error = [ K1mean-K1min, K1max-K1mean]
##axs[0].errorbar(range(31),K1mean, yerr = asymmetric_error)
##axs[0].plot(K1_pha_obs, 'r-*', label = 'observation')
##K1mean=np.average(K1_amp, axis=1)
##K1max = np.max(K1_amp,axis=1)
##K1min = np.min(K1_amp,axis=1)
##asymmetric_error = [ K1mean-K1min, K1max-K1mean]
##axs[1].errorbar(range(31),K1mean, yerr = asymmetric_error)
##axs[1].plot(K1_amp_obs, 'r-*', label = 'observation')



