#
# Script for DETIDING
# by Salish Sea MEOPAR
# https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis/raw/tip/compare_tides/Analysis8Components.ipynb
#
print("Setting the environment...")
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
print("..Done!")
###
#### pathname for data - all of the tide runs are stored in this directory
####path = '/data/nsoontie/MEOPAR/SalishSea/results/tides/'
path = '/work/ag15419/tmp/diaharm_ana/'
print("Work dir is",path )
###
####the run we want to analyze
####runname = 'corr15'
runname = 'simt_ctrl0'
print("Runname is",runname )
###
####joining the two string together
#name = path +runname +'/'
name = path
###
print("Name:",name)
###
stations = 'Gibraltar'
print("Stations list: ",stations)
###
numsta=len(stations)
print("Stations num:",numsta )
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
print("M2freq",M2freq,"radians per hour")
print("S2freq",S2freq,"radians per hour")

#### initial phase calculation
M2ft = 1.036581
M2uvt = 0.009395

S2ft = 1.000000
S2uvt = 0.000000

N2ft = 1.037550
N2uvt = 0.00464

K1ft = 0.882850
K1uvt =  

print("Initial phases")
print("M2 ft and uvt:",M2ft,M2uvt)
print("S2 ft and uvt:",S2ft,S2uvt)
###
def double(x, M2amp, M2pha, S2amp, S2pha):
    return (M2amp*np.cos(M2freq*x-M2pha*np.pi/180.)+
            S2amp*np.cos(S2freq*x-S2pha*np.pi/180.))

###
#fig, ax = plt.subplots(1,1,figsize=(12,5))
#for stn in (0,4,14,23):
#    print(stations[stn])
#    fT1 = NC.Dataset(name+stations[stn]+'.nc','r')
#    time = fT1.variables["time_counter"][:]/3600. # want hours not seconds
#    ssh = fT1.variables["sossheig"][:,0,0]
#    ax.plot(time,ssh)
print(stations)
fT1 = NC.Dataset(name+stations+'.nc','r')
time = fT1.variables["time_counter"][:]/3600. # want hours not seconds
ssh = fT1.variables["sossheig"][:,0,0]
#print(ssh)
#fig = plt.figure(figsize=(9,5.5))
#ax = fig.add_subplot(111)
#ax.plot(time,ssh)
#plt.savefig('ssh.jpg')
###

#print =(ssh.shape)
print(ssh)
print(time)
print(ssh.shape[0])
###
M2_amp=[]; M2_pha=[]; S2_amp=[]; S2_pha=[]
###
#ts = 150
te = ssh.shape[0]
#print(ts)
print(te)
#print("time-step:",ts)
###
###for stn in range(numsta):
#ncfilename=name+stations+'.nc'
#print(name)
#print(stations)
#print(ncfilename)

#fT1 = NC.Dataset(name+stations+'.nc','r')
#print(fT1.data_model)
#time = fT1.variables["time_counter"][:]/3600. # want hours not seconds
#ssh = fT1.variables["sossheig"][:,0,0]
###
#print(ssh)
ts=0
print(time[ts:te])
fitted, cov = curve_fit(double,time[ts:te],ssh[ts:te])
print(fitted[:])
if fitted[0] < 0:
    fitted[0] = -fitted[0]
    fitted[1] = fitted[1]+180
M2_amp.append(fitted[0]*M2ft)
pha = fitted[1]+M2uvt
if pha > 360:
   pha=pha-360
elif pha < 0:
   pha = pha+360
###if stn == 6:
###print pha
M2_pha.append(pha)
###
if fitted[2] < 0:
   fitted[2] = -fitted[2]
   fitted[3] = fitted[3]+180
S2_amp.append(fitted[2]*S2ft)
pha= fitted[3]+S2uvt
if pha > 360:
   pha=pha-360
S2_pha.append(pha)
###
###
# PLOT
plt.figure(figsize=(20,12))
###
#plt.subplot(1,2,1)
plt.plot(np.array(M2_amp)*1.03, '-bo', label = 'model')
####plt.plot(M2_amp_obs, 'r-o', label = 'observation')
plt.title('M2 Amplitude')
#plt.legend( loc='upper left' )
###
###
#plt.subplot(1,2,2)
#plt.plot(np.array(S2_amp)*0.99, '-bo', label = 'model')
#plt.title('S2 Amplitude')

plt.savefig('M2_S2_A.jpg')
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
