#!/bin/bash
#
# ACG 23/06/2019
# Script for hourly TS extraction  
# Ini file: p_extr.ini 
#
#set -u
set -e
#set -x 
################### PREPROC ###################################

# Source ini file
  source ts_extr.ini

# Read and check infos (work dir, file names, archive dir, etc.)

  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     cp ${SRC_DIR}/ts_extr.ini ${ANA_WORKDIR}
     
     # Clean workdir
     echo "WARNING: I am going to remove all files in $ANA_WORKDIR ..."
     sleep 10
     for TO_BE_RM in $( ls $ANA_WORKDIR ); do
         #rm $ANA_WORKDIR/$TO_BE_RM
         echo $TO_BE_RM
     done

  else
     echo "ERROR: WORKDIR $ANA_WORKDIR NOT FOUND!!"
     exit
  fi

  # Input file check
  
  # Num check
  if [[ ${#ANA_INPATHS[@]} -ne ${#ANA_INFILES_TPL[@]} ]]; then
     echo "ERROR: Check Inputs array, something is missing!"
     exit
  else
     INSET_NUM=${#ANA_INPATHS[@]}
     echo "The num of input file sets is $INSET_NUM"
  fi

  # File check

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    if [[ -d ${ANA_INPATHS[${IDX_IN}]} ]]; then

      IDX_DATE=$ANA_STARTDATE
      
      while [[ $IDX_DATE -le $ANA_ENDDATE ]]; do

        echo "Date: $IDX_DATE"
        ANA_INFILE=$( echo ${ANA_INFILES_TPL[${IDX_IN}]} | sed -e s/%YYYYMMDD%/$IDX_DATE/g  )

        FOUND_NUM=0
        if [[ -e ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE ]]; then
        #if [[ -e ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE ]]; then
           FOUND_NUM=$(( $FOUND_NUM + 1 ))
           echo "Found infile: $ANA_INFILE"
           ln -sf ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE .
           #ln -sf ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE .
        else
           echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/"
        fi
        #echo "I have found $FOUND_NUM files in ${ANA_INPATHS[$IDX_IN]} dataset.."

      IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d ) 
      done    
 
    else 
      echo "ERROR: Input dir ${ANA_INPATHS[${IDX_IN}]} NOT FOUND!!"
      exit
    fi

    IDX_IN=$(( $IDX_IN + 1 ))
    done

# Read ana type from ini file and set the environment
  
  if [[ $TS_FLAG == 1 ]]; then
    echo "TS analisys for points listed in $TS_COOFILE is required!"
    module load $TS_MODULE
    echo "Enviroment.."
  fi


####################### ANALISYS ##############################
# +-----------------------+
# | POINT                 |
# +-----------------------+ 

# Read point coo ( and set job array )
    
     echo "Reading STZ coordinates.."

     declare -a COO_TS_LAT
     declare -a COO_TS_LON
     declare -a COO_TS_NAME

     echo "Reading STZ coordinates.."
     COO_IDX=0
     while read COOFILE_LINE; do
       if [[ ${COOFILE_LINE:0:1} != "#" ]]; then
          COO_TS_LAT[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 1 -d";" )
          COO_TS_LON[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 2 -d";" )
          COO_TS_NAME[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 3 -d";" )
          echo "STZ NAME: ${COO_TS_NAME[${COO_IDX}]} "
          echo "LAT/LON: ${COO_TS_LAT[${COO_IDX}]}/${COO_TS_LON[${COO_IDX}]}"
          if [[ $OBS_FLAG == 1 ]]; then
            OBS_PATH[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 4 -d";" )
            OBS_NAME[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 5 -d";" )
            echo "OBS_FILES: ${OBS_PATH[$COO_IDX]}/${OBS_NAME[$COO_IDX]}"
          elif [[ $OBS_FLAG == 2 ]]; then
            OBS_PATH[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 4 -d";" )
            OBS_NAME[${COO_IDX}]=$( echo $COOFILE_LINE | cut -f 5 -d";" )
            echo "OBS_FILES: ${OBS_PATH[$COO_IDX]}/${OBS_NAME[$COO_IDX]}"
          fi
       COO_IDX=$(( $COO_IDX + 1 ))
       fi 
     done < $TS_COOFILE

     echo "COO NUM = ${#COO_TS_NAME}"

# Extract values from netCDF

 if [[ $MOD_FLAG == 1 ]]; then

  echo "------ TS Extraction from model outputs ------"

  TS_OUTFILE_PRE="ts"
  TS_OUTFILE_TPL="${TS_OUTFILE_PRE}_${ANA_STARTDATE}_${ANA_ENDDATE}_%STZ%_%FIELD%_%IDX%.txt"

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    echo "I am workin on ${ANA_INTAG[$IDX_IN]} dataset..."

    for VAR in ${TS_FIELDS[@]}; do

     echo "Extracting field: $VAR .. "

     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )
     
     TS_IDX=0
     while [[ $TS_IDX -lt $COO_IDX ]]; do

      NC_TS_OUTFILE_TPL="%STZ%_mod_${ANA_INTAG[$IDX_IN]}.nc"
      NC_TS_OUTFILE=$( echo "$NC_TS_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )

      echo "# TS extraction from ${ANA_INTAG[$IDX_IN]}" 
      echo "# STZ: ${COO_TS_NAME[$TS_IDX]}" 
      echo "# VAR: ${VAR}" 
      echo "# OUTFILE = $NC_TS_OUTFILE" 
 
      # Extraction
      IDX_NC=0
      for TSEXT in $EXT_INS; do

        # TEMPORARY (NEW)
        LON_TOCUT=${COO_TS_LON[$TS_IDX]}
        LAT_TOCUT=${COO_TS_LAT[$TS_IDX]}
        LON_INT=$( echo $LON_TOCUT | cut -f 1 -d"." )
        LAT_INT=$( echo $LAT_TOCUT | cut -f 1 -d"." )
        cdo sellonlatbox,$(( ${LON_INT} - 1 )),$(( ${LON_INT} + 1 )),$(( ${LAT_INT} - 1 )),$(( ${LAT_INT} + 1 )) $TSEXT red_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo setctomiss,0.0000 red_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc 


        # OLD ----
        # Point selection With land/sea mask (Nan)
        #cdo setctomiss,0.0000 $TSEXT miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc #Old ok
        #cdo -remapdis,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}_${IDX_NC}.nc #Old ok
        #------------
         
        # TEMPORARY (NEW)
        cdo setmisstonn miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc nnmiss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo -remapnn,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} nnmiss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}_${IDX_NC}.nc
        IDX_NC=$(( ${IDX_NC} + 1 ))
      done

  

        # Mean comp
        cdo mergetime tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}_*.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo timmean tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo sub tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        cdo mergetime zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}*.nc 30_mod.nc
        if [[ $OBS_FLAG == "1" ]] || [[ $OBS_FLAG == "0" ]] ; then
           cdo intntime,2 30_mod.nc 00_mod.nc 
           cdo seltime,00:00:00,01:00:00,02:00:00,03:00:00,04:00:00,05:00:00,06:00:00,07:00:00,08:00:00,09:00:00,10:00:00,11:00:00,12:00:00,13:00:00,14:00:00,15:00:00,16:00:00,17:00:00,18:00:00,19:00:00,20:00:00,21:00:00,22:00:00,23:00:00,24:00:00 00_mod.nc $NC_TS_OUTFILE
           rm -v 00_mod.nc 
           rm -v 30_mod.nc
        elif [[ $OBS_FLAG == "2" ]] ; then
           mv 30_mod.nc $NC_TS_OUTFILE
           rm -v 30_mod.nc
        fi
        #rm -v 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        rm -v miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        rm -v red_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc # new
        rm -v nnmiss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc # new
        rm -v mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        rm -v tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}*.nc
        rm -v zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

     TS_IDX=$(( $TS_IDX + 1 ))
     done

    done

   IDX_IN=$(( $IDX_IN + 1 ))
   done

 fi


  if [[ $OBS_FLAG == 1 ]]; then
  
     echo "Obs extraction.."

# Extract values from netCDF TG obs

    echo "------ OBS Extraction ------"

    for O_VAR in ${OBS_VAR[@]}; do 

      echo "Obs Var: $O_VAR"
         
        STZ_IDX=0
        while [[ $STZ_IDX -lt $COO_IDX ]]; do
          O_STZ=${COO_TS_NAME[$STZ_IDX]}
          O_PATH=${OBS_PATH[$STZ_IDX]}
          O_NAME=${OBS_NAME[$STZ_IDX]}
       
         if [[ $O_PATH != "NaN" ]] && [[ $O_NAME != "NaN" ]] ; then

          NCOBS_OUTFILE_TPL="%STZ%_obs.nc"
          NCOBS_OUTFILE=$( echo "$NCOBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
          MODNC_OUTFILE_TPL='diff_modnc_%STZ%.nc'
          NCOBS_OUTFILE=$( echo "$NCOBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
       
          echo "# Obs data " 
          echo "# TS in ${O_STZ} " 
          echo "# VAR: ${O_VAR}" 
          echo "# DT: ${ANA_STARTDATE:0:6} - ${ANA_ENDDATE:0:6}" 

           ID_DATE=$ANA_STARTDATE
           while [[ ${ID_DATE:0:6} -le ${ANA_ENDDATE:0:6} ]]; do
             OBS_FILE=$( echo "${O_PATH}/${O_NAME}" | sed -e "s/%YYYY%/${ID_DATE:0:4}/g" -e "s/%YYYYMM%/${ID_DATE:0:6}/g" )


             for OBS_FOUND in $( ls $OBS_FILE ); do
              if [[ -f $OBS_FOUND ]]; then
                 cp $OBS_FOUND ncobs_tmp_${ID_DATE}.nc
                 cdo expr,"sossheig=SLEV*1" ncobs_tmp_${ID_DATE}.nc sossobs_tmp_${ID_DATE}.nc
              else
                echo "NOT found OBS file: $OBS_FILE "
              fi
             done
            

           ID_DATE=$( date -d "${ID_DATE} 1 month" +%Y%m%d )
           done

           # Merge tmp files, select field and change var name
             
             cdo mergetime sossobs_tmp_*.nc merged_obs.nc
             #cdo expr,"sossheig=SLEV*1" merged_obs.nc tosub_obs.nc
             cdo timmean merged_obs.nc obs_mean.nc
             cdo sub merged_obs.nc obs_mean.nc $NCOBS_OUTFILE
             #cdo expr,"sossheig=SLEV*1" merged_obs.nc $NCOBS_OUTFILE
             # Time interpolation (for comparison with ISPRA obs and mod )
             ###cdo intntime,2 00_obs.nc $NCOBS_OUTFILE
             rm -v ncobs_tmp_*.nc
             rm -v merged_obs.nc
             #rm -v tosub_obs.nc
             rm -v obs_mean.nc
             ###rm -v 00_obs.nc
     
             # Compute diffs between mod and nc obs, stat values, spectra and edf
             IDX_IN=0
             while [[ $IDX_IN -lt $INSET_NUM ]]; do
                  # Model inputs
                  NC_TS_OUTFILE_TPL="%STZ%_mod_${ANA_INTAG[$IDX_IN]}.nc"
                  NC_TS_OUTFILE=$( echo "$NC_TS_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
  
                  # Diffs
                  DIFF_MODNN_TPL="diff_modnc_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  DIFF_MODNN=$( echo "$DIFF_MODNN_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )

                  cdo sub $NC_TS_OUTFILE $NCOBS_OUTFILE $DIFF_MODNN

                  # Mean and bias
                  MOD_MEAN_TPL="mean_mod_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  MOD_MEAN=$( echo "$MOD_MEAN_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
                  cdo timmean  $NC_TS_OUTFILE $MOD_MEAN

                  NCOBS_MEAN_TPL="mean_ncobs_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  NCOBS_MEAN=$( echo "$NCOBS_MEAN_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
                  cdo timmean  $NCOBS_OUTFILE $NCOBS_MEAN

                  NCOBSMOD_BIAS_TPL="bias_modncobs_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  NCOBSMOD_BIAS=$( echo "$NCOBSMOD_BIAS_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
                  cdo sub  $MOD_MEAN $NCOBS_MEAN $NCOBSMOD_BIAS
                 

                  # RMSE

                  # SPECTRUM

                  #MOD_SPT_TPL="spect_mod_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  #MOD_SPT=$( echo "$MOD_SPT_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
                  #cdo fourier,1  $NC_TS_OUTFILE $MOD_SPT

                  #NCOBS_SPT_TPL="mean_ncobs_%STZ%_${ANA_INTAG[$IDX_IN]}.nc"
                  #NCOBS_SPT=$( echo "$NCOBS_SPT_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
                  #cdo fourier  $NCOBS_OUTFILE $NCOBS_SPT
                  
                  # EDF       
  

             IDX_IN=$(( $IDX_IN + 1 ))
             done
          fi
         STZ_IDX=$(( $STZ_IDX + 1 ))
         done         
        
    done

  elif [[ $OBS_FLAG == 2 ]]; then

    echo "------ ISPRA OBS Extraction ------"

    for O_VAR in ${OBS_VAR[@]}; do

      echo "Obs Var: $O_VAR"

         STZ_IDX=0
         while [[ $STZ_IDX -lt $COO_IDX ]]; do
         O_STZ=${COO_TS_NAME[$STZ_IDX]}
         O_PATH=${OBS_PATH[$STZ_IDX]}
         O_NAME=${OBS_NAME[$STZ_IDX]}

         if [[ $O_PATH != "NaN" ]] && [[ $O_NAME != "NaN" ]] ; then

          OBS_OUTFILE_PRE="obs"
          OBS_OUTFILE_TPL="${OBS_OUTFILE_PRE}_%STZ%_${ANA_STARTDATE}_${ANA_ENDDATE}.csv"

          OBS_OUTFILE=$( echo "$OBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
          #echo "# Obs data - TS in ${O_STZ} - VAR: ${O_VAR} - DT: ${ANA_STARTDATE:0:8} - ${ANA_ENDDATE:0:8}" > $OBS_OUTFILE
          #echo "# TS in ${O_STZ} " >> $OBS_OUTFILE
          #echo "# VAR: ${O_VAR}" >> $OBS_OUTFILE
          #echo "# DT: ${ANA_STARTDATE:0:8} - ${ANA_ENDDATE:0:8}" >> $OBS_OUTFILE
          #echo "# " >> $OBS_OUTFILE
          echo "idNum;station;year;month;day;hour;minu;sec;value" > $OBS_OUTFILE
           ID_LINE=1
           ID_DATE=${ANA_STARTDATE:0:8}
           ISPRA_ENDDATE=$( date -d "${ANA_ENDDATE:0:8} 1 day" +%Y%m%d )
           while [[ ${ID_DATE:0:8} -le ${ISPRA_ENDDATE} ]]; do
             echo "DATE: ${ID_DATE:0:8} "
             OBS_FILE=$( echo "${O_PATH}/${O_NAME}" | sed -e "s/%YYYYMMDD%/${ID_DATE:0:8}/g" )

             for OBS_FOUND in $( ls $OBS_FILE ); do
              if [[ -f $OBS_FOUND ]]; then
                while read ISPRA_LINE ; do
                   if [[ ${ISPRA_LINE:0:1} != "G" ]]; then
                      ISPRA_DATA=${ISPRA_LINE:0:8}
                      ISPRA_ORA=$( echo $ISPRA_LINE | cut -f 2 -d";" )
                      ISPRA_VAL=$( echo $ISPRA_LINE | cut -f 3 -d";" | sed -e "s/","/"."/g" )
                      if [[ ${ISPRA_ORA:3:2} == "30"  ]]; then
                         echo $( echo "$ID_LINE;${O_STZ};20${ISPRA_DATA:6:2};${ISPRA_DATA:3:2};${ISPRA_DATA:0:2};${ISPRA_ORA:0:2};${ISPRA_ORA:3:2};00;${ISPRA_VAL}" )  >> $OBS_OUTFILE
                         ID_LINE=$(( $ID_LINE + 1 ))
                      fi
                   fi
                done < $OBS_FOUND
                #cat $OBS_FILE >> $OBS_OUTFILE
              else
                echo "NOT found OBS file: $OBS_FILE "
              fi
             done

           ID_DATE=$( date -d "${ID_DATE} 1 day" +%Y%m%d )
           done
          fi
         STZ_IDX=$(( $STZ_IDX + 1 ))
         done

    done
                                                                                                                                                       
  fi

###################### POSTPROC ###########################

# Output check

# Archive


