#!/bin/bash
#
#########################################################################
# This script reads (and writes in a csv file) from exp files
# found in /work/%user%/exp dirs the following infos  
# 1) exp name
# 2) start date, duration
# 3) active keys 
# 4) namelist and iofile name and path  
# 5) computational resources infos
#
# Then it extracts from nemo-med-simulations 
# 1) the global wall time (which includes pre and post processing phases)
# from log.txt files
# 2) the NEMO stand-alone wall-time from bhist_P? files computing a mean
# on daily values (since daily run are performed)
# 3) Output file.nc dimensions
#
# Written by AC Goglio on 14/06/2019
#######################################################################
#######################################
# USERS SETTINGS:
MY_WORK_EXP="/work/ag15419/exp/eas5"
MY_EXP_SRC="/users/home/ag15419/src_dev/eas5_2"

CSV_OUTFILE="/users/home/ag15419/tides_pp/exp_pp/foundExp_$( date -u +%Y%m%d%H%M%S).csv"
echo "CSV_OUTFILE=$CSV_OUTFILE"
echo "# This csv file lists the nemo-med exp infos"
echo "# for exps found in $MY_WORK_EXP"
echo "# "
echo "# Exp_name;start_date;duration;restart_flag;NEMO_keys;namelist_path(existence_flag);iofile_path;iofile_fields;comp_resources;mean_Drun_walltime;outfile_dims;" > $CSV_OUTFILE
######################################
######################################
#
# Found EXP
EXP=$( ls ${MY_WORK_EXP} )
echo "In ${MY_WORK_EXP} there are the following exp:"
echo $EXP

echo "####"
EXP_NUM=0
for EXP_NAME in $EXP ; do 
    DESCR_FILE="${MY_WORK_EXP}/${EXP_NAME}/exp-descriptor.sh"
    if [[ -e ${DESCR_FILE} ]]; then
       EXP_NUM=$(( $EXP_NUM + 1 ))
       echo "#$EXP_NUM) $EXP_NAME"
       echo "Found ${EXP_NAME} descriptor file..."
       # Start date from descriptor file
       START_DATE=$( grep "timing_start_time" ${DESCR_FILE} | cut -f 2 -d"'" )
       echo "START_DATE=$START_DATE"
       # Duration must be computed from outfile 1d/T availability
       DURATION_DESCRIPTOR=$( grep "timing_hours" ${DESCR_FILE} | cut -f 2 -d"'" | cut -f 1 -d"#" )
       DURATION=$( ls -t ${MY_WORK_EXP}/${EXP_NAME}/output/*/*.nc | grep "_1d_" | grep "_T" | wc -l )
       echo "DURATION=$DURATION [h]"
       # Restart flag from comparison between desciptor timing hours and outfile availability
       DURATION_DESCR_DAYS=$(( ${DURATION_DESCRIPTOR} / 24 ))
       #RESTAR_FLAG=$( grep "timing_start_from_restart" | cut -f 2 -d"=" )
       if [[ ${DURATION_DESCR_DAYS} -gt $DURATION ]]; then
          RESTAR_FLAG=1
       else
          RESTAR_FLAG=0
          if [[ ${DURATION_DESCR_DAYS} -lt $DURATION ]]; then
             echo "WARNING: check run exit status..."
          fi
       fi 
       echo "RESTAR_FLAG=$RESTAR_FLAG"
       # NEMO keys (PROB: is there a copy of the file with the list of the active keys in the work dir? )
       KEYS="??"
       echo "KEYS=$KEYS"
       # Namelist name from descriptor file
       NAMELIST_NAME=$( grep "NEMO_NL" ${DESCR_FILE} | cut -f 2 -d"=" )
       NAMELIST_PATH="${MY_EXP_SRC}/pack/med_prod/phys/${NAMELIST_NAME}"
       echo "NAMELIST=$NAMELIST_NAME"
       if [[ -e $NAMELIST_PATH ]];then
          NAMELIST_FLAG=1 
       else
          NAMELIST_FLAG=0
       fi
       echo "NAMELIST_FLAG=$NAMELIST_FLAG"
       # iodef (field selection)
       IODEF_NAME=$( grep "EXP_LF" ${DESCR_FILE} | cut -f 2 -d"=" )
       echo "IODEF_NAME=$IODEF_NAME"
       # work copy of iodef file
       MY_TMP_IODEF=${MY_WORK_EXP}/${EXP_NAME}/tmp/iodef.xml
       if [[ -e ${MY_TMP_IODEF} ]]; then
          echo -n $( sed -n '/TRUE/,/FALSE/p' ${MY_TMP_IODEF} | grep -v "FALSE" | grep "long_name" | cut -f 4 -d"\"" ) 
       fi
       # NEMO time-step
       NEMO_TIMESTEP=$( grep "NEMOTimestep" ${DESCR_FILE} | cut -f 2 -d"=" )
       echo "NEMO_TIMESTEP=$NEMO_TIMESTEP"
       

       # Computational resources
       RESOURCE_PPN=$( grep "nemo_n_mpi_proc" ${DESCR_FILE} | cut -f 2 -d"=" | grep -v "#" )
       echo "RESOURCE_PPN=$RESOURCE_PPN"
    else
       echo "NOT Found ${DESCR_FILE} descriptor file..." 
    fi
done
echo "####"
#####################################
# TIME AND DIMS
#####################################
# Pack log files path and name

LOF_EXP=$EXP
LOF_PATH=${MY_WORK_EXP}
LOF_NAME="log.txt"

BHIST_PRENAME="bhist_P"


START_DATE=( 20160101 20160101 20150101 20160101 20160101 )

# Loop on exps
DATE_ID=0
for EXP_ID in ${LOF_EXP}; do
    EXP_START_DATE=${START_DATE[$DATE_ID]}
    echo "################################"
    echo -e "EXP = ${EXP_ID}\n"
    
    LOF=${LOF_PATH}/${EXP_ID}/${LOF_NAME} 
    ##echo "Log file: $LOF"
    BHISTS="${LOF_PATH}/${EXP_ID}/output/${EXP_START_DATE:0:6}/$BHIST_PRENAME*"
    #echo "bhist files: $BHISTS"
    #
    #######################################################
    # Tot DT
    START_TIME=$( grep "_init" $LOF | cut -f 10 -d" " )
    END_TIME=$( grep "End of all" $LOF | cut -f 9 -d" " )
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
    ##echo "$START_TIME --> $END_TIME" 
    #echo "H1 M1 S1 = $H1 $M1 $S1 --> H2 M2 S2 = $H2 $M2 $S2"

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
