#!/bin/bash
#BSUB -J map_extr          # Name of the job.
#BSUB -o /work/ag15419/job_scratch/maps_%J.out  # Appends std output to file %J.out.
#BSUB -e /work/ag15419/job_scratch/maps_%J.err  # Appends std error to file %J.err.
#BSUB -cwd "/work/ag15419/job_scratch/%J/"
#BSUB -q serial_24h
#BSUB -n 1    # Number of CPUs
#
#
# ACG 04/11/2019
# Script for map TS extraction 
# Ini file: map_extr.ini 
#
#set -u
set -e
#set -x 
################### ENV SETTINGS ##############################
module load CDO
#cdo cat /work/ag15419/tmp/currents/nc_extr/name_tmp_*.nc /work/ag15419/tmp/currents/nc_extr/all_tmp.nc
cdo timmean /work/ag15419/tmp/currents/nc_extr/all_tmp.nc  /work/ag15419/tmp/currents/nc_extr/map3D_yr_allv_sossheig_20160101_20181231_mod_Tides8.nc
cdo timmean /work/ag15419/tmp/currents/nc_extr_ctrl/all_tmp.nc  /work/ag15419/tmp/currents/nc_extr_ctrl/map3D_yr_allv_sossheig_20160101_20181231_mod_Control.nc
rm -v /work/ag15419/tmp/currents/nc_extr/all_tmp.nc
rm -v /work/ag15419/tmp/currents/nc_extr_ctrl/all_tmp.nc
###################### POSTPROC ###########################

# Output check

# Archive


