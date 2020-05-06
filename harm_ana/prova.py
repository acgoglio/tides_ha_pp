# imports
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
import plotly
from plotly import graph_objects as go # for bar plot
#
#######################################################
path = '/work/ag15419/tmp/h_ana_miss/'
iniend_dates= '01/07/2017 - 31/12/2017'
dates_lab='20170701_20171231'
name = path
########################################################
# Global plots

plt.figure(figsize=(12,6))
plt.rc('font', size=12)
    #
#plt.figure(
#   # data=[
#  go.Bar(
#            name=GLOB_A_mod[0][1],
#            x=GLOB_A_mod[:][0],
#            y=GLOB_A_mod[:][1],
#            #offsetgroup=0,
#        ),
#  go.Bar(
#            name=GLOB_A_mod[0][2],
#            x=GLOB_A_mod[:][0],
#            y=GLOB_A_mod[:][2],
#            #offsetgroup=1,
#            #base=GLOB_A_mod[:][1],
#        ),
#  go.Bar(
#            name=GLOB_A_obs[0][1],
#            x=GLOB_A_obs[:][0],
#            y=GLOB_A_obs[:][1],
#            #\offsetgroup=0,
#        ),
#  go.Bar(
#            name=GLOB_A_obs[0][2],
#            x=GLOB_A_obs[:][0],
#            y=GLOB_A_obs[:][2],
#            #offsetgroup=1,
#            #base=GLOB_A_obs[:][1],
#        )
#    #],
#    #layout=go.Layout(
#    #    title='Tidal Amplitudes MODEL Vs OBS'+'('+iniend_dates+')',
#    #    yaxis_title="Tidal Components Amplitudes [cm]"
#    )

data = {
    "original":[15, 23, 32, 10, 23],
    "model_1": [4,   8, 18,  6,  0],
    "model_2": [11, 18, 18,  0,  20],
    "labels": [
        "feature",
        "question",
        "bug",
        "documentation",
        "maintenance"
    ]
}


fig = go.Figure(
    data=[
        plt.Bar(
            name="Original",
            x=data["labels"],
            y=data["original"],
            offsetgroup=0,
        ),
        plt.Bar(
            name="Model 1",
            x=data["labels"],
            y=data["model_1"],
            offsetgroup=1,
        ),
        plt.Bar(
            name="Model 2",
            x=data["labels"],
            y=data["model_2"],
            offsetgroup=1,
            base=data["model_1"],
        )
    ],
    layout=plt.Layout(
        title="Issue Types - Original and Models",
        yaxis_title="Number of Issues"
    )
)
#offl.plot({'data': data, 'layout': layout}, filename='GLOBAL_A_'+dates_lab+'.jpg', auto_open=False)

#fig.show()
#fig.write_image('GLOBAL_A_'+dates_lab+'.jpg')
plt.savefig(path+'GLOBAL_A_'+dates_lab+'.jpg')
#plt.clf()
