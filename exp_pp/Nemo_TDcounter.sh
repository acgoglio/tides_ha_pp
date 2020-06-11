#!/bin/bash
#
###############################################
# This script extracts from nemo-med-simulations 
# 1) the global wall time (which includes pre and post processing phases)
# from log.txt files
# 2) the NEMO stand-alone wall-time from bhist_P? files computing a mean
# on daily values (since daily run are performed)
# 3) Output file.nc dimensions
#
#
# Written by AC Goglio on 12/06/2019
##############################################


#############################################
# USERS SETTINGS:

# Exp names ( must be a string with gaps between names )
#LOF_EXP="simu_ctrl0 simt_nt_2015 simt_nt_16 simt_nt_15 simu_ctrl0_ldread simt_nt75_16" 
LOF_EXP="simu_ctrl0_NewIntScheme simu_ctrl0_NewIntScheme1 simu_ctrl0_NewIntScheme2 simu_ctrl0_NewIntScheme3"
LOF_EXP="simt_tra3 simt_tra simt_tra2 simt_ctrl0_NObdy"
#LOF_EXP="simu_tides8_6"
# Exp start date ( must be an array of values formatted as yyyymmdd, one element per exp ) 
#START_DATE=( 20160101 20150101 20160101 20150101 20160101 20160101 )
START_DATE=( 20181008 20150101 20170128 20180928 )
START_DATE=( 20150101 20150101 20150101 20150101 )
#START_DATE=( 20150101 )

######################################

#Pack log files path and name
LOF_PATH="/work/ag15419/exp/eas5"
#LOF_PATH="/work/ec04916/exp/eas5/"
LOF_NAME="log.txt"

BHIST_PRENAME="bhist_P"

# Loop on exps
DATE_ID=0
for EXP_ID in ${LOF_EXP}; do
    EXP_START_DATE=${START_DATE[$DATE_ID]}
    echo "################################"
    echo -e "EXP = ${EXP_ID}\n"
    
    LOF=${LOF_PATH}/${EXP_ID}/${LOF_NAME} 
    ##echo "Log file: $LOF"
    BHISTS="${LOF_PATH}/${EXP_ID}/log/${EXP_START_DATE:0:6}/$BHIST_PRENAME*" # output (old system) or log (new)
    #echo "bhist files: $BHISTS"
    #
    #######################################################
    # Tot DT
    START_TIME=$( grep "_init" $LOF | tail -n 1 | cut -f 10 -d" " )
    END_TIME=$( grep "End of all" $LOF | tail -n 1 | cut -f 9 -d" " )
    # Hours
    H1=$( echo $START_TIME | cut -f 1 -d\: )
    H2=$( echo $END_TIME | cut -f 1 -d\: )
    # Mins
    M1=$( echo $START_TIME | cut -f 2 -d\: )
    M2=$( echo $END_TIME | cut -f 2 -d\: )
    # Secs
    S1=$( echo $START_TIME | cut -f 3 -d\: )
    S2=$( echo $END_TIME | cut -f 3 -d\: )
    ##echo -e "\n\n\n+--------------------------+"    
    ##echo "TOT time"
    echo "$START_TIME --> $END_TIME" 
    echo "H1 M1 S1 = $H1 $M1 $S1 --> H2 M2 S2 = $H2 $M2 $S2"

    # DELTA computing
     for T_ID in H1 H2 M1 M2 S1 S2 ; do
       #echo -e "\n${T_ID}=${!T_ID}"
       if [[ "${!T_ID:0:1}" == "0" ]]; then
          ##echo "First char of $T_ID = 0, I am going to remove it.."
          declare "${T_ID}"=${!T_ID:1:1}
          ##echo -e "\t...Done!"
       #else
          #new_${T_ID}=${!T_ID}
       fi
       #echo -e "${T_ID}=${!T_ID}"
     done
    
    DELTA_h=$(( ${H2} - ${H1} ))
    #echo "DH=$DELTA_h"
    DELTA_m=$(( ${M2} - ${M1} ))
    #echo "DM=$DELTA_m"
    DELTA_s=$(( ${S2} - ${S1} ))
    #echo "DS=$DELTA_s"
    
    TOT_sec=$(( $DELTA_h * 3600 + DELTA_m * 60 + $DELTA_s ))
    ##echo "********************"
    echo "TOT sec = $TOT_sec s"
    ##echo "********************"

    # NEMO wall-time
    ##echo -e "\n\n\n+--------------------------+"    
    ##echo "NEMO wall-time"
    # NEMO log file
    ##OCE_LOF="${LOF_PATH}/${EXP_ID}/output/${EXP_START_DATE:0:6}/timing.output_*00"
    ##for OCE_LOF_ID in $( ls $OCE_LOF ); do
    ##    grep -A 4 "Total timing (sum)" $OCE_LOF_ID | tail -n 2 | cut -f 1
    ##done

    # BHIST time
    COUNT=0
    TOT_RUN_TIME=0
    for BHIST_FOUND in $( ls $BHISTS ); do
        RUN_TIME=$( grep -A 2 "Summary"  $BHIST_FOUND | grep -v "\-\-" | grep -v "Summary" | grep -v "PEND" | awk '{print $3}' )
        ##grep -A 2 "Summary"  $BHIST_FOUND | grep -v "\-\-" | grep -v "Summary" | grep -v "PEND" | awk '{print $3}' 
        ##echo "RUN_TIME=$RUN_TIME"
        TOT_RUN_TIME=$(( $TOT_RUN_TIME + $RUN_TIME ))
        COUNT=$(( $COUNT + 1 ))
    done
    ##echo $COUNT
    ##echo $TOT_RUN_TIME
    if [[ $COUNT != 0 ]]; then
       MEAN_RUNTIME=$(( $TOT_RUN_TIME / $COUNT ))
       echo "Nemo daily wall-time=$MEAN_RUNTIME s"
    fi
    ################################
    # Output dimensions
    NEMO_OUTFILES="${LOF_PATH}/${EXP_ID}/output/${EXP_START_DATE:0:6}/*${EXP_START_DATE:0:6}*.nc"
    NEMO_FIRSTDAY="${LOF_PATH}/${EXP_ID}/output/${EXP_START_DATE:0:6}/*${EXP_START_DATE}*.nc"
    NETCDF_NUM=$( ls -l ${NEMO_OUTFILES} | wc -l )
     if [[ $NETCDF_NUM != 0 ]]; then
        echo "# of nc files = $NETCDF_NUM"
        TOT_DIM=$( du -hc $NEMO_OUTFILES | grep "total" )
        echo "Tot dim: $TOT_DIM"
        FIRSTDAY_DIM=$( du -hc $NEMO_FIRSTDAY | grep "total" )
        echo "First day dim: $FIRSTDAY_DIM"
     else 
        echo "Empty output dir: WHY????"
     fi
DATE_ID=$(( $DATE_ID + 1 ))
done
