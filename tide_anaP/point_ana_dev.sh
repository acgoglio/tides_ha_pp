#!/bin/bash
#
# ACG 23/06/2019
# Script for field comparisons  
# Ini file: tide_ana.ini 
#
#set -u
set -e
#set -x 
################### PREPROC ###################################

# Source ini file
  source tide_ana.ini

# Read and check infos (work dir, file names, archive dir, etc.)

if [[ $EXTRACTION_FLAG == 1 ]] || [[ $DMMM_FLAG == 1 ]] ; then

  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     cp ${SRC_DIR}/tide_ana.ini ${ANA_WORKDIR}/tide_ana.ini_$(date +%Y%m%d_%H%M%S)
     
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

 if [[ $MYDIAG_FLAG == 0 ]]; then
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
 fi

else
  echo "Non-extraction run.." 
  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     cp ${SRC_DIR}/tide_ana.ini ${ANA_WORKDIR}
  else
     echo "ERROR: WORKDIR $ANA_WORKDIR NOT FOUND!!"
     exit

  fi
fi

# Read ana type from ini file and set the environment
  
  if [[ $DMMM_FLAG == 1 ]]; then
    echo "Min, mean and max analisys over the whole domain is required!"
    module load $DMMM_MODULE
  fi

  if [[ $TS_FLAG == 1 ]]; then
    echo "TS analisys for points listed in $TS_COOFILE is required!"
    module load $TS_MODULE
  fi

  if [[ $MYDIAG_FLAG != 0 ]]; then
    echo "Mydiag analysis environment!"
    module load $MYDIAG_MODULE
  fi

 

####################### ANALISYS ##############################

 if [[ $DMMM_FLAG == 1 ]]; then
 
# +-----------------------+
# | AREA                  |
# +-----------------------+


  ##########################
# DMMM Extraction and plotting
 
  echo "------ Extraction ------"
   
  DMMM_OUTFILE_PRE="dmmm"
  DMMM_OUTFILE_TPL="${DMMM_OUTFILE_PRE}_${ANA_STARTDATE}_${ANA_ENDDATE}_%FIELD%_%IDX%.txt"

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    echo "I am workin on ${ANA_INTAG[$IDX_IN]} dataset..."

    for VAR in ${DMMM_FIELDS[@]}; do 
 
      echo "Extracting field: $VAR .. "
 
      EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )
      DMMM_OUTFILE=$( echo "$DMMM_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g"  )

      echo "# DMMM extraction form ${ANA_INTAG[$IDX_IN]}" > $DMMM_OUTFILE 
      echo "# " >> $DMMM_OUTFILE
      echo "# -1 :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name" >> $DMMM_OUTFILE

      # Extraction
      for DMMEXT in $EXT_INS; do 
        # Standard:
        echo "cdo infon $DMMEXT | grep -v "Date" | grep "$VAR"  >> $DMMM_OUTFILE"
        cdo infon $DMMEXT | grep -v "Date" | grep "$VAR"  >> $DMMM_OUTFILE
        #
        # Quasi-Atlantic_box cut:
        #echo "cdo -infon -sellonlatbox,-18.125,-6.000,30.1875,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_1" 
        #cdo -infon -sellonlatbox,-18.125,-6.000,30.1875,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_1
        #echo "cdo -infon -sellonlatbox,-6.000,1.000,45.000,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_2" 
        #cdo -infon -sellonlatbox,-6.000,1.000,42.500,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_2
        # Quasi Med cut:
        #echo "cdo -infon -sellonlatbox,1.000,36.29167,30.1875,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_1" 
        #cdo -infon -sellonlatbox,1.000,36.29167,30.1875,45.97917 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_1
        #echo "cdo -infon -sellonlatbox,-5.990,0.000,30.1875,42.000 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_2" 
        #cdo -infon -sellonlatbox,-5.990,0.000,30.1875,42.000 $DMMEXT | grep -v "Date" | grep "$VAR"  >> ${DMMM_OUTFILE}_2
        #
        #rm -v $DMMEXT
      done

     done

    IDX_IN=$(( $IDX_IN + 1 ))
    done 

    #########################
    echo "------ Plot ------"
    # DMMM Plots

    IDX_PVAR=0

    for PVAR in ${DMMM_FIELDS[@]}; do

      echo "Plotting field: $PVAR .. "
   
      #TXT_TABS=$( echo "$DMMM_OUTFILE_TPL" | sed -e "s/%IDX%/"*"/g" -e "s/%FIELD%/${VAR}/g"  )

      PLOT_VAR=$PVAR
      PLOT_UNITS=${DMMM_UDM[$IDX_PVAR]}

      PLOT_SDATE=$ANA_STARTDATE
      PLOT_EDATE=$ANA_ENDDATE

      DMMM_PLOT=$( echo "$DMMM_PLOTFILE_TPL" | sed -e "s/%DATES%/"${PLOT_SDATE}_${PLOT_EDATE}"/g" -e "s/%FIELD%/${PVAR}/g"  )
      echo "Outplot: $DMMM_PLOT"

      if [[ $INSET_NUM -eq 2 ]]; then
         
        echo "Comparison plot with 2 datasets..."

        TXT_TAB_0=$( echo "$DMMM_OUTFILE_TPL" | sed -e "s/%IDX%/0/g" -e "s/%FIELD%/${PVAR}/g" )
        TXT_TAB_1=$( echo "$DMMM_OUTFILE_TPL" | sed -e "s/%IDX%/1/g" -e "s/%FIELD%/${PVAR}/g" )

        GNUPLOT_TMP="tmp.gpl_${PVAR}"


#        cat << EOF > ${GNUPLOT_TMP}  

        echo "#" > ${GNUPLOT_TMP}
        echo "set term jpeg size 1000,600 giant " >> ${GNUPLOT_TMP}
        echo "set output \"$DMMM_PLOT\" " >> ${GNUPLOT_TMP}
        echo "set title \"Max, mean and min on MED-MFC Domain ; VAR: $PVAR  ; DT: $ANA_STARTDATE - $ANA_ENDDATE\" " >> ${GNUPLOT_TMP}
        echo "set xlabel \"Time\" " >> ${GNUPLOT_TMP}
        echo "set xdata time " >> ${GNUPLOT_TMP}
        echo "set timefmt \"%Y-%m-%d %H:%M:%S\" " >> ${GNUPLOT_TMP}
        echo "set xrange [\"${PLOT_SDATE:0:4}-${PLOT_SDATE:4:2}-${PLOT_SDATE:6:2} 00:00:00\":\"${PLOT_EDATE:0:4}-${PLOT_EDATE:4:2}-${PLOT_EDATE:6:2} 23:30:00\"] " >> ${GNUPLOT_TMP}
        echo "set format x \"%d/%m/%Y\" " >> ${GNUPLOT_TMP}
        echo "set ylabel \"$PLOT_VAR ${PLOT_UNITS}\" " >> ${GNUPLOT_TMP}
        echo "set grid " >> ${GNUPLOT_TMP}
        echo "plot '$TXT_TAB_0' using 3:9 with line lw 2 title \"Min ${ANA_INTAG[0]}\" , '$TXT_TAB_0' using 3:10 with line lw 2 title \"Mean ${ANA_INTAG[0]}\", '$TXT_TAB_0' using 3:11 with line lw 2 title \"Max ${ANA_INTAG[0]}\", '$TXT_TAB_1' using 3:9 with line lw 2 title \"Min ${ANA_INTAG[1]}\" , '$TXT_TAB_1' using 3:10 with line lw 2 title \"Mean ${ANA_INTAG[1]}\", '$TXT_TAB_1' using 3:11 with line lw 2 title \"Max ${ANA_INTAG[1]}\" " >> ${GNUPLOT_TMP}
        #echo "plot '$TXT_TAB_0' using 3:10 with line lw 2 title \"Mean ${ANA_INTAG[0]}\",'$TXT_TAB_1' using 3:10 with line lw 2 title \"Mean ${ANA_INTAG[1]}\" " >> ${GNUPLOT_TMP}

#EOF 

        gnuplot < $GNUPLOT_TMP
     fi
  
   IDX_PVAR=$(( $IDX_PVAR + 1 ))
   done
 fi
# +-----------------------+
# | POINT                 |
# +-----------------------+ 

# Read point coo ( and set job array )
 if [[ $TS_FLAG == 1 ]]; then

     declare -a COO_TS_LAT
     declare -a COO_TS_LON
     declare -a COO_TS_NAME

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
if [[ $EXTRACTION_FLAG == 1 ]]; then

    echo "------ TS Extraction ------"

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

      TS_OUTFILE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
 
      echo "PROVA: # TS extraction from ${ANA_INTAG[$IDX_IN]} in $TS_OUTFILE"
      echo "# TS extraction from ${ANA_INTAG[$IDX_IN]}" > $TS_OUTFILE
      echo "# " >> $TS_OUTFILE
      echo "# STZ: ${COO_TS_NAME[$TS_IDX]}" >> $TS_OUTFILE
      echo "# VAR: ${VAR}" >> $TS_OUTFILE
      echo "# " >> $TS_OUTFILE
      echo "# -1 :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name" >> $TS_OUTFILE

      # Extraction
      IDX_NC=0
      for TSEXT in $EXT_INS; do
     
        # Time interpolation (for comparison with obs)
        cdo intntime,2 $TSEXT 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc 
        # With land/sea mask (Nan)
        cdo setctomiss,0.0000 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo -remapdis,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}_${IDX_NC}.nc
        IDX_NC=$(( ${IDX_NC} + 1 ))
      done
  
        # Mean comp
        cdo mergetime tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}_*.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        if [[ $DETR_FLAG == 1 ]]; then
           cdo timmean tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
           cdo sub tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
           rm -v mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
           rm -v tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}*.nc
        else 
           mv tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        fi

        rm -v 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        rm -v miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        cdo -infon zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc  | grep "$VAR" >> $TS_OUTFILE 
        rm -v zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc


     TS_IDX=$(( $TS_IDX + 1 ))
     done

    done

   IDX_IN=$(( $IDX_IN + 1 ))
   done


  if [[ $OBS_FLAG == 1 ]]; then

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

          OBS_OUTFILE_PRE="obs_ts"
          OBS_OUTFILE_TPL="${OBS_OUTFILE_PRE}_${ANA_STARTDATE}_${ANA_ENDDATE}_%FIELD%_%STZ%.txt"

          OBS_OUTFILE=$( echo "$OBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
          echo "# Obs data " > $OBS_OUTFILE
          echo "# TS in ${O_STZ} " >> $OBS_OUTFILE
          echo "# VAR: ${O_VAR}" >> $OBS_OUTFILE
          echo "# DT: ${ANA_STARTDATE:0:6} - ${ANA_ENDDATE:0:6}" >> $OBS_OUTFILE
          echo "# " >> $OBS_OUTFILE      

           ID_DATE=$ANA_STARTDATE
           while [[ ${ID_DATE:0:6} -le ${ANA_ENDDATE:0:6} ]]; do
             echo "DATE: ${ID_DATE:0:6} ${ID_DATE} ${ANA_ENDDATE:0:6}"
             OBS_FILE=$( echo "${O_PATH}/${O_NAME}" | sed -e "s/%YYYY%/${ID_DATE:0:4}/g" -e "s/%YYYYMM%/${ID_DATE:0:6}/g" )


             # with zero sub
             for OBS_FOUND in $( ls $OBS_FILE ); do
              if [[ -f $OBS_FOUND ]]; then
                 if [[ $DETR_FLAG == 1 ]]; then
                    cdo timmean $OBS_FILE obs_mean.nc
                    cdo sub $OBS_FILE obs_mean.nc obs_zero.nc
                    cdo infon obs_zero.nc | grep "${O_VAR} " >> $OBS_OUTFILE
                    rm -v obs_mean.nc
                    rm -v obs_zero.nc
                 else
                    mv $OBS_FILE obs_zero.nc
                    cdo infon obs_zero.nc | grep "${O_VAR} " >> $OBS_OUTFILE
                    rm -v obs_zero.nc
                 fi
              else
                echo "NOT found OBS file: $OBS_FILE "
              fi
             done
             # without zero 
             #cdo infon $OBS_FILE | grep "${O_VAR} " >> $OBS_OUTFILE

           ID_DATE=$( date -d "${ID_DATE} 1 month" +%Y%m%d )
           done
          fi
         STZ_IDX=$(( $STZ_IDX + 1 ))
         done         
        
    done

  # Extraction of ISPRA TG data from txt 
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

          OBS_OUTFILE_PRE="obs_ts"
          OBS_OUTFILE_TPL="${OBS_OUTFILE_PRE}_${ANA_STARTDATE}_${ANA_ENDDATE}_%FIELD%_%STZ%.txt"

          OBS_OUTFILE=$( echo "$OBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${O_STZ}/g" )
          echo "# Obs data " > $OBS_OUTFILE
          echo "# TS in ${O_STZ} " >> $OBS_OUTFILE
          echo "# VAR: ${O_VAR}" >> $OBS_OUTFILE
          echo "# DT: ${ANA_STARTDATE:0:8} - ${ANA_ENDDATE:0:8}" >> $OBS_OUTFILE
          echo "# " >> $OBS_OUTFILE

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
                      ISPRA_VAL=$( echo $ISPRA_LINE | cut -f 3 -d";" )

                      #ISPRA_DATASTR="20${ISPRA_DATA:6:2}${ISPRA_DATA:3:2}${ISPRA_DATA:0:2}"
                      #echo "ISPRA_DATASTR=$ISPRA_DATASTR DATE=${ID_DATE}"
                      #if [[ ${ID_DATE} == ${ISPRA_DATASTR} ]] && [[ $DIFF_FLAG == 0 ]]; then
                      #   echo " $ID_LINE : 20${ISPRA_DATA:6:2}-${ISPRA_DATA:3:2}-${ISPRA_DATA:0:2} ${ISPRA_ORA}:00       0        1       0 :                 ${ISPRA_VAL}    : sossheig "  >> $OBS_OUTFILE
                      #elif [[ ${ID_DATE} == ${ISPRA_DATASTR} ]] && [[ $DIFF_FLAG == 1 ]]; then
                      #   echo " $ID_LINE : 20${ISPRA_DATA:6:2}-${ISPRA_DATA:3:2}-${ISPRA_DATA:0:2} ${ISPRA_ORA}:00       0        1       0 :                 ${ISPRA_VAL}    : sossheig " | grep " 00:00 " >> $OBS_OUTFILE
                      #fi

                      if [[ $DIFF_FLAG == 0 ]]; then
                         echo " $ID_LINE : 20${ISPRA_DATA:6:2}-${ISPRA_DATA:3:2}-${ISPRA_DATA:0:2} ${ISPRA_ORA}:00       0        1       0 :                 ${ISPRA_VAL}    : sossheig "  >> $OBS_OUTFILE
                      elif [[ $DIFF_FLAG == 1 ]] && [[ ${ISPRA_ORA:3:2} == "00" ]] ; then
                         echo " $ID_LINE : 20${ISPRA_DATA:6:2}-${ISPRA_DATA:3:2}-${ISPRA_DATA:0:2} ${ISPRA_ORA}:00       0        1       0 :                 ${ISPRA_VAL}    : sossheig " >> $OBS_OUTFILE
                      fi
                         
                   ID_LINE=$(( $ID_LINE + 1 ))
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
# End of EXTRACTION_FLAG 
fi


# Plot ts

  echo "------ Plot TS ------"

    echo "I am workin on # ${INSET_NUM} dataset(s)..."

    ID_V=0
    for VAR in ${TS_FIELDS[@]}; do
      O_VAR=${OBS_VAR[${ID_V}]}
      ID_V=$(( $ID_V + 1 ))  

      IDX_UDMVAR=0
      echo "Plot field: $VAR .. "

      PLOT_UNITS=${DMMM_UDM[$IDX_PVAR]}
      PLOT_SDATE=$ANA_STARTDATE
      PLOT_EDATE=$ANA_ENDDATE


      TS_IDX=0
      while [[ $TS_IDX -lt $COO_IDX ]]; do
      echo "STZ: ${COO_TS_NAME[$TS_IDX]}"

        TSPLOT_INFILE_0=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/0/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
        TSPLOT_INFILE_1=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/1/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
        TSPLOT_INFILE_2=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/2/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
 
        if [[ $DIFF_FLAG == 1 ]] && [[ $OBS_FLAG == 0 ]]; then
          # Case to be checked...
          echo "DIFF active"
          INFILE_PASTE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/all/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
            paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_1 > $INFILE_PASTE
            echo "INFILE_PASTE=$INFILE_PASTE"         
        fi

        if [[ $OBS_FLAG != 0 ]]; then
           TSPLOT_INFILE_OBS=$( echo "$OBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
           if [[ $DIFF_FLAG == 1 ]] && [[ $OBS_FLAG == 2 ]]; then
              echo "DIFF FLAG ON!"

              INFILE_PASTE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/allobs/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )

              if [[ ${MOD1_FLAG} == 1 ]]; then
                 # Case to be checked..
                 #paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_OBS > $INFILE_PASTE
                 while read OLINE ; do
                   if [[ ${OLINE:0:1} != "#" ]]; then
                      echo $OLINE
                      ODATA=$( echo $OLINE | cut -f 1 -d\ )
                      echo "DATA: $ODATA"
                      OTIME=$( echo $OLINE | cut -f 2 -d\ )
                      echo "TIME: $OTIME"
                      OVAL=$( echo $OLINE | cut -f 3 -d\ )
                      echo "OVAL: $OVAL"
                      echo "$( grep -c "${ODATA} ${OTIME}" $TSPLOT_INFILE_0 )"
                      EXIST_DATATIME=$( echo "$( grep -c "${ODATA} ${OTIME}" $TSPLOT_INFILE_0 )" )
                      echo "EXIST_DATATIME=$EXIST_DATATIME"
                      if [[ $EXIST_DATATIME == 1 ]]; then
                         MVAL=$( grep "$ODATA_TIME" $TSPLOT_INFILE_0 )
                         echo "$MVAL $OVAL " >> $INFILE_PASTE
                      fi
                   fi
                 done < $TSPLOT_INFILE_OBS
              else
                 paste $TSPLOT_INFILE_1 $TSPLOT_INFILE_2 > paste_temp.txt
                 grep -v "#" $TSPLOT_INFILE_OBS > ${TSPLOT_INFILE_OBS}_clean
                 CLEAN_LINE=$( wc -l ${TSPLOT_INFILE_OBS}_clean )
                 if [[ ${CLEAN_LINE:0:1} != 0 ]]; then
                   awk '{print $3,$4,$9 }' ${TSPLOT_INFILE_OBS}_clean > ${TSPLOT_INFILE_OBS} 
                   rm -v ${TSPLOT_INFILE_OBS}_clean
                   #paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_1 $TSPLOT_INFILE_OBS > $INFILE_PASTE
                   while read OLINE ; do
                     if [[ "${OLINE:0:1}" != "#" ]]; then
                        #ODATA_TIME=$( echo $OLINE | cut -f 1 -d\  )
                        ODATA=$( echo $OLINE | cut -f 1 -d\ )
                        OTIME=$( echo $OLINE | cut -f 2 -d\ )
                        OVAL=$( echo $OLINE | cut -f 3 -d\ )
                        EXIST_DATATIME=$( echo "$( grep -c "${ODATA} ${OTIME}" paste_temp.txt )" )
                        if [[ $EXIST_DATATIME == 1 ]]; then
                            MVALS=$( grep "${ODATA} ${OTIME}" paste_temp.txt )
                            echo "$MVALS $OVAL " >> $INFILE_PASTE
                        fi
                     fi
                   done < $TSPLOT_INFILE_OBS
                 fi
                 rm -v paste_temp.txt
              fi

           elif [[ $DIFF_FLAG == 1 ]] && [[ $OBS_FLAG == 1 ]]; then
              echo "DIFF FLAG ON! (Obs from netCDF)"
              INFILE_PASTE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/allobs/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
              if [[ ${MOD1_FLAG} == 1 ]]; then
                 echo "To be written.."
              else
                 paste $TSPLOT_INFILE_1 $TSPLOT_INFILE_2 > paste_temp.txt
                 while read OLINE ; do
                   if [[ "${OLINE:0:1}" != "#" ]]; then
                      ODATA=$( echo $OLINE | cut -f 3 -d\ )
                      OTIME=$( echo $OLINE | cut -f 4 -d\ )
                      OVAL=$( echo $OLINE | cut -f 5 -d: )
                      EXIST_DATATIME=$( echo "$( grep -c "${ODATA} ${OTIME}" paste_temp.txt )" )
                      if [[ $EXIST_DATATIME == 1 ]]; then
                          MVALS=$( grep "${ODATA} ${OTIME}" paste_temp.txt )
                          echo "$MVALS $OVAL " >> $INFILE_PASTE
                      fi
                   fi
                 done < $TSPLOT_INFILE_OBS

                 #rm -v paste_temp.txt

              fi
           fi
        fi

        TS_PLOT=$( echo "$TS_PLOT_TPL" | sed -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" -e "s/%DATES%/"${PLOT_SDATE}_${PLOT_EDATE}"/g" )
       
        if [[ $DIFF_FLAG == 0 ]]; then 

          # Plot mod ts 
          GNUPLOT_TS_TMP="ts_tmp.gpl"    

          echo "#" > ${GNUPLOT_TS_TMP}
          echo "set term jpeg size 1500,600 giant" >> ${GNUPLOT_TS_TMP}
          echo "set output \"$TS_PLOT\" " >> ${GNUPLOT_TS_TMP}
          if [[ $OBS_FLAG == 2 ]]; then
             echo "stats '$TSPLOT_INFILE_OBS' using 9 nooutput" >> ${GNUPLOT_TS_TMP}
          fi
          echo "set title \"Time Series (P: ${COO_TS_NAME[$TS_IDX]}  VAR: $VAR  DT: $ANA_STARTDATE - $ANA_ENDDATE ) \"  " >> ${GNUPLOT_TS_TMP}
          # 1 model..
          #echo "stats '$TSPLOT_INFILE_0' using 9 prefix \"Mod\", '$TSPLOT_INFILE_OBS' using 9 prefix \"Obs\" " >> ${GNUPLOT_TS_TMP}
          #
          echo "set xlabel \"Time\" " >> ${GNUPLOT_TS_TMP}
          echo "set xdata time " >> ${GNUPLOT_TS_TMP}
          echo "set timefmt \"%Y-%m-%d %H:%M:%S\" " >> ${GNUPLOT_TS_TMP}
          echo "set xrange [\"${PLOT_SDATE:0:4}-${PLOT_SDATE:4:2}-${PLOT_SDATE:6:2} 00:00:00\":\"${PLOT_EDATE:0:4}-${PLOT_EDATE:4:2}-${PLOT_EDATE:6:2} 23:30:00\"] " >> ${GNUPLOT_TS_TMP}
          echo "set format x \"%d/%m/%Y\" " >> ${GNUPLOT_TS_TMP}
          echo "set ylabel \"$PLOT_VAR ${PLOT_UNITS}\" " >> ${GNUPLOT_TS_TMP}
          echo "set grid " >> ${GNUPLOT_TS_TMP}
          if [[ $OBS_FLAG == 0 ]]; then
             echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
          elif [[ $OBS_FLAG == 1 ]]; then
             echo "Obs plot TG netCDF"

             #3 mod
             if [[ ${#ANA_INTAG[@]} == 3 ]]; then
               echo "3 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\", '$TSPLOT_INFILE_OBS' using 3:9 with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}

             #2 mod
             elif [[ ${#ANA_INTAG[@]} == 2 ]]; then
               echo "2 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_OBS' using 3:9 with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}
           
             # 1 model..
             elif [[ ${#ANA_INTAG[@]} == 1 ]]; then
               echo "1 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_OBS' using 3:9 with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}
             fi
             #echo "plot '$TSPLOT_INFILE_OBS' using 3:9 with line lw 2 title \"Obs\" " >> ${GNUPLOT_TS_TMP}

          elif [[ $OBS_FLAG == 2 ]]; then
             echo "Obs plot ISPRA TG txt"
             if [[ $DETR_FLAG == 1 ]]; then
                #3 mod
                if [[ ${#ANA_INTAG[@]} == 3 ]]; then
                  echo "3 mod"
                  echo "plot '$TSPLOT_INFILE_OBS' using 3:(\$9/100 - STATS_mean/100 )  with line lw 2 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\" " >> ${GNUPLOT_TS_TMP}

                #2 mod
                elif [[ ${#ANA_INTAG[@]} == 2 ]]; then
                  echo "2 mod"
                  echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100 - STATS_mean/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}

                # 1 model..
                elif [[ ${#ANA_INTAG[@]} == 1 ]]; then
                  echo "1 mod"
                  echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100-STATS_mean/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}
                fi

             # NO Detrending
             else
                #3 mod
                if [[ ${#ANA_INTAG[@]} == 3 ]]; then
                  echo "3 mod"
                  echo "plot '$TSPLOT_INFILE_OBS' using 3:(\$9/100)  with line lw 2 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\" " >> ${GNUPLOT_TS_TMP}

                #2 mod
                elif [[ ${#ANA_INTAG[@]} == 2 ]]; then
                  echo "2 mod"
                  echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}

                # 1 model..
                elif [[ ${#ANA_INTAG[@]} == 1 ]]; then
                  echo "1 mod"
                  echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}
                fi

             fi
          fi

        else 
         if [[ $DETR_FLAG == 0 ]]; then 
          # TS + Diff (only for detrended case..)


          # Plot mod ts 
          GNUPLOT_TS_TMP="ts_tmp.gpl"

          echo "#" > ${GNUPLOT_TS_TMP}
          echo "set term jpeg size 1500,900 giant" >> ${GNUPLOT_TS_TMP}
          echo "set output \"$TS_PLOT\" " >> ${GNUPLOT_TS_TMP}
          if [[ $OBS_FLAG == 2 ]]; then
             echo "stats '$TSPLOT_INFILE_OBS' using 3 nooutput" >> ${GNUPLOT_TS_TMP}
             echo "stats '$INFILE_PASTE' using (\$9 - (\$23/100 - STATS_mean/100)) name 'STATS1' nooutput" >> ${GNUPLOT_TS_TMP}
             echo "stats '$INFILE_PASTE' using (\$20 - (\$23/100 - STATS_mean/100)) name 'STATS2' nooutput" >> ${GNUPLOT_TS_TMP}
          elif [[ $OBS_FLAG == 1 ]]; then
             echo "stats '$TSPLOT_INFILE_OBS' using 9 nooutput" >> ${GNUPLOT_TS_TMP}
             echo "stats '$INFILE_PASTE' using (\$9 - (\$23 - STATS_mean)) name 'STATS1' nooutput" >> ${GNUPLOT_TS_TMP}
             echo "stats '$INFILE_PASTE' using (\$20 - (\$23 - STATS_mean)) name 'STATS2' nooutput" >> ${GNUPLOT_TS_TMP}
          fi

          echo "set multiplot layout 2,1 title \"Time Series (P: ${COO_TS_NAME[$TS_IDX]}  VAR: $VAR  DT: $ANA_STARTDATE - $ANA_ENDDATE )\" " >> ${GNUPLOT_TS_TMP}
    
          echo "set xlabel \"Time\" " >> ${GNUPLOT_TS_TMP}
          echo "set xdata time " >> ${GNUPLOT_TS_TMP}
          echo "set timefmt \"%Y-%m-%d %H:%M:%S\" " >> ${GNUPLOT_TS_TMP}
          echo "set xrange [\"${PLOT_SDATE:0:4}-${PLOT_SDATE:4:2}-${PLOT_SDATE:6:2} 00:00:00\":\"${PLOT_EDATE:0:4}-${PLOT_EDATE:4:2}-${PLOT_EDATE:6:2} 23:30:00\"] " >> ${GNUPLOT_TS_TMP}
          echo "set format x \"%d/%m/%Y\" " >> ${GNUPLOT_TS_TMP}
          echo "set ylabel \"$PLOT_VAR ${PLOT_UNITS}\" " >> ${GNUPLOT_TS_TMP}
          echo "set grid " >> ${GNUPLOT_TS_TMP}
          echo "set xzeroaxis lw 2" >> ${GNUPLOT_TS_TMP}
          if [[ $OBS_FLAG == 0 ]]; then
             echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
             echo "plot '$INFILE_PASTE' using 3:(\$9-\$20) with line lw 2 title \"${ANA_INTAG[0]}-${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
          elif [[ $OBS_FLAG == 2 ]]; then
             echo "Obs plot"
             echo "plot '$TSPLOT_INFILE_OBS' using 1:(\$3/100 - STATS_mean/100 )  with points pt 7 ps 0.5 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\" " >> ${GNUPLOT_TS_TMP}
             # Colors and order for NO bdy case..
             #echo "plot '$TSPLOT_INFILE_OBS' using 1:(\$3/100 - STATS_mean/100 )  with points pt 7 ps 0.5 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"green\" title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
            echo "plot '$INFILE_PASTE' using 3:(\$20 - (\$23/100 - STATS_mean/100) ) with points pt 1 ps 0.5 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS\",'$INFILE_PASTE' using 3:(\$9 - (\$23/100 - STATS_mean/100) ) with points pt 1 ps 0.5 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS\", STATS1_mean lw 2 lc rgb \"black\" title gprintf(\"${ANA_INTAG[1]}-OBS mean %g [${PLOT_UNITS}] \", STATS1_mean) , STATS1_lo_quartile lw 2 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS 1st quartile\",STATS1_up_quartile lw 2 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS 3rd quartile\",STATS2_mean lw 2 lc rgb \"blue\" title gprintf( \"${ANA_INTAG[2]}-OBS mean = %g [${PLOT_UNITS}] \", STATS2_mean ), STATS2_lo_quartile lw 2 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS 1st quartile\",STATS2_up_quartile lw 2 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS 3rd quartile\" " >> ${GNUPLOT_TS_TMP}
          #echo " plot \"$INFILE_PASTE\" u (hist((\$9 - (\$23/100 - STATS_mean/100)),width)):(1.0)/STATS1_records smooth freq w boxes lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS err distribution\", plot \"$INFILE_PASTE\" u (hist((\$20 - (\$23/100 - STATS_mean/100)),width)):(1.0)/STATS2_records smooth freq w boxes lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS err distribution\" " >> ${GNUPLOT_TS_TMP}
          elif [[ $OBS_FLAG == 1 ]]; then
             echo "Obs plot"
             echo "plot '$TSPLOT_INFILE_OBS' using 3:(\$9 - STATS_mean )  with points pt 7 ps 0.5 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\" " >> ${GNUPLOT_TS_TMP}
             # Colors and order for NO bdy case..
             #echo "plot '$TSPLOT_INFILE_OBS' using 3:(\$9 - STATS_mean ) with points pt 7 ps 0.5 lt rgb \"blue\" title \"Obs\" , '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"green\" title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
            echo "plot '$INFILE_PASTE' using 3:(\$20 - (\$23 - STATS_mean) ) with points pt 1 ps 1 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS\",'$INFILE_PASTE' using 3:(\$9 - (\$23 - STATS_mean) ) with points pt 1 ps 1 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS\", STATS1_mean lw 2 lc rgb \"black\" title gprintf(\"${ANA_INTAG[1]}-OBS mean %g [${PLOT_UNITS}] \", STATS1_mean) , STATS1_lo_quartile lw 2 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS 1st quartile\",STATS1_up_quartile lw 2 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS 3rd quartile\",STATS2_mean lw 2 lc rgb \"blue\" title gprintf( \"${ANA_INTAG[2]}-OBS mean = %g [${PLOT_UNITS}] \", STATS2_mean ), STATS2_lo_quartile lw 2 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS 1st quartile\",STATS2_up_quartile lw 2 lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS 3rd quartile\" " >> ${GNUPLOT_TS_TMP}
          #echo " plot \"$INFILE_PASTE\" u (hist((\$9 - (\$23/100 - STATS_mean/100)),width)):(1.0)/STATS1_records smooth freq w boxes lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS err distribution\", plot \"$INFILE_PASTE\" u (hist((\$20 - (\$23/100 - STATS_mean/100)),width)):(1.0)/STATS2_records smooth freq w boxes lc rgb \"red\" title \"${ANA_INTAG[2]}-OBS err distribution\" " >> ${GNUPLOT_TS_TMP}
          fi
         else
          echo "Computing the diff between non detrended time-series does not make sense!I"
         fi
        fi

        # Plot
        gnuplot < $GNUPLOT_TS_TMP || echo "Prob with this plot..why?!"
        mv -v $GNUPLOT_TS_TMP ${GNUPLOT_TS_TMP}_${COO_TS_NAME[$TS_IDX]}
        #rm -v $GNUPLOT_TS_TMP

      TS_IDX=$(( $TS_IDX + 1 ))
      done

     IDX_UDMVAR=$(( $IDX_UDMVAR + 1 ))
     done

 fi


# +-----------------------+
# | MYDIAG DIFF           |
# +-----------------------+ 

 if [[ $MYDIAG_FLAG == 1 ]]; then
    echo "Mydiag analysis.."

    IDX_MYDIAG=0
    for MYDIAG_FIELD in ${MYDIAG_SHORT_NAMES[@]}; do
      #if [[ $IDX_MYDIAG != 2 ]] && [[ $IDX_MYDIAG != 5 ]] && [[ $IDX_MYDIAG != 8 ]] && [[ $IDX_MYDIAG != 11 ]] && [[ $IDX_MYDIAG != 14 ]] && [[ $IDX_MYDIAG != 25 ]] && [[ $IDX_MYDIAG != 33 ]] && [[ $IDX_MYDIAG != 35 ]]; then
        echo "$IDX_MYDIAG ) I am working on ${MYDIAG_SHORT_NAMES[${IDX_MYDIAG}]}.."

        # Infiles
        MYDIAG_FILE_1="${MYDIAG_PATH_1}/${MYDIAG_FILES[${IDX_MYDIAG}]}.nc"
        MYDIAG_FILE_2="${MYDIAG_PATH_2}/${MYDIAG_FILES[${IDX_MYDIAG}]}.nc"
         
        #     
        if [[ -f $MYDIAG_FILE_1 ]] && [[ -f $MYDIAG_FILE_2 ]] ; then

          # Select dates 
          cdo seldate,${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-${ANA_STARTDATE:6:2}T00:00:00,${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}T23:30:00 $MYDIAG_FILE_1 seldate1.nc
          cdo seldate,${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-${ANA_STARTDATE:6:2}T00:00:00,${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}T23:30:00 $MYDIAG_FILE_2 seldate2.nc

          # Compute the difference 
          MYDIAG_DIFFOUT="Diff_${MYDIAG_FIELD}.nc"
          cdo sub seldate1.nc seldate2.nc $MYDIAG_DIFFOUT 

          # Extract TS
          MYDIAG_PLOTINFILE_1="plotinfile_1_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_2="plotinfile_2_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_DIFF="plotinfile_diff_${MYDIAG_FIELD}.txt"

          # The following two line should be RM 
          #cdo outputtab,date,time,value,name $MYDIAG_FILE_1 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_1 
          #cdo outputtab,date,time,value,name $MYDIAG_FILE_2 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_2
          # 
          cdo outputtab,date,time,value,name seldate1.nc | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_1 
          cdo outputtab,date,time,value,name seldate2.nc | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_2

          cdo outputtab,date,time,value,name $MYDIAG_DIFFOUT | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_DIFF
         

          # Plot TS and Diff
          PLOT_SDATE=$ANA_STARTDATE
          PLOT_EDATE=$ANA_ENDDATE

          MYDIAG_PLOTUDM=${MYDIAG_UDM[${IDX_MYDIAG}]} 

          MYDIAG_PLOT=$( echo "$MYDIAG_PLOT_TPL" | sed -e "s/%FIELD%/${MYDIAG_FILES[${IDX_MYDIAG}]}/g"  -e "s/%DATES%/"${PLOT_SDATE}_${PLOT_EDATE}"/g" )
          echo "MYDIAG_PLOT_TPL $MYDIAG_PLOT"


          GNUPLOT_MYDIAG_TMP="mydiag_tmp.gpl_${MYDIAG_FILES[${IDX_MYDIAG}]}"
          #GNUPLOT_MYDIAG_TMPTXT="mydiag_tmptxt.gpl"


          ## Txt gpl file
          #echo "#" > ${GNUPLOT_MYDIAG_TMPTXT}
          #echo "set term jpeg giant size 1800,900 font \"Times,45\"" >> ${GNUPLOT_MYDIAG_TMPTXT}
          #echo "set output \"$MYDIAG_PLOT\" " >> ${GNUPLOT_MYDIAG_TMPTXT}

          #echo "stats '$MYDIAG_PLOTINFILE_1' using 3 name 'STATS1'" >> ${GNUPLOT_MYDIAG_TMPTXT}
          #echo "stats '$MYDIAG_PLOTINFILE_2' using 3 name 'STATS2'" >> ${GNUPLOT_MYDIAG_TMPTXT}
          #echo "stats '$MYDIAG_PLOTINFILE_DIFF' using 3 name 'STATSD'" >> ${GNUPLOT_MYDIAG_TMPTXT}

          # Plot gpl file
          echo "#" > ${GNUPLOT_MYDIAG_TMP}
          echo "set term jpeg size 1200,600 giant" >> ${GNUPLOT_MYDIAG_TMP} #1300,800
          echo "set output \"$MYDIAG_PLOT\" " >> ${GNUPLOT_MYDIAG_TMP}

          echo "stats '$MYDIAG_PLOTINFILE_1' using 3 name 'STATS1' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_2' using 3 name 'STATS2' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_DIFF' using 3 name 'STATSD' nooutput" >> ${GNUPLOT_MYDIAG_TMP}

          echo "set multiplot layout 2,1 title \"Time Series ( VAR: ${MYDIAG_LONG_NAMES[${IDX_MYDIAG}]}  DT: $ANA_STARTDATE - $ANA_ENDDATE )\" " >> ${GNUPLOT_MYDIAG_TMP}

          echo "set key opaque" >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xlabel \"Time\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xdata time " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set timefmt \"%Y-%m-%d %H:%M:%S\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xrange [\"${PLOT_SDATE:0:4}-${PLOT_SDATE:4:2}-${PLOT_SDATE:6:2} 00:00:00\":\"${PLOT_EDATE:0:4}-${PLOT_EDATE:4:2}-${PLOT_EDATE:6:2} 23:30:00\"] " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set format x \"%d/%m/%Y\" " >> ${GNUPLOT_MYDIAG_TMP}
          #echo "set ylabel \"${MYDIAG_FIELD} ${MYDIAG_PLOTUDM}\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set grid " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set key left top" >> ${GNUPLOT_MYDIAG_TMP} # right bottom
          echo "set xzeroaxis lt 2 lc rgb \"black\" lw 3" >> ${GNUPLOT_MYDIAG_TMP}

          # All:
          #echo "plot '$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 3 lt rgb \"${MYDIAG_COLOR1}\" title gprintf(\"${MYDIAG_INTAG_1} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS1_mean), '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 3 lt rgb \"${MYDIAG_COLOR2}\" title gprintf(\"${MYDIAG_INTAG_2} (AVG = %g ${MYDIAG_PLOTUDM} )\", STATS2_mean) " >> ${GNUPLOT_MYDIAG_TMP}
          # All no avg:
          echo "set ylabel \"${MYDIAG_LONG_NAMES[${IDX_MYDIAG}]} ${MYDIAG_PLOTUDM}\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "plot '$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 3 lt rgb \"${MYDIAG_COLOR1}\" title \"${MYDIAG_INTAG_1}\", '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 3 lt rgb \"${MYDIAG_COLOR2}\" title \"${MYDIAG_INTAG_2}" >> ${GNUPLOT_MYDIAG_TMP} 
          # SSH:
          #echo "plot STATS1_mean lw 2 lc rgb \"cyan\" title gprintf( \"AVG ${MYDIAG_INTAG_1}=%g ${MYDIAG_PLOTUDM}\", STATS1_mean ), STATS2_mean lw 2 lc rgb \"gray\" title gprintf( \"AVG ${MYDIAG_INTAG_2}=%g ${MYDIAG_PLOTUDM}\", STATS2_mean ),'$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 2 lt rgb \"blue\" title gprintf(\"${MYDIAG_INTAG_1}\", STATS1_mean), '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 2 lt rgb \"black\" title gprintf(\"${MYDIAG_INTAG_2}\", STATS2_mean)" >> ${GNUPLOT_MYDIAG_TMP}
          # The following 2 lines should be REMOVED (written for compaison with Cucco 2016, transprts at Messina strait comparison)
          #echo "set yrange [\"-0.5\":\"0.6\"] " >> ${GNUPLOT_MYDIAG_TMP}
          #echo "plot '$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 2 lt rgb \"black\" title gprintf(\"${MYDIAG_INTAG_1} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS1_mean), '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 2 lt rgb \"red\" title gprintf(\"${MYDIAG_INTAG_2} (AVG = %g ${MYDIAG_PLOTUDM} )\", STATS2_mean) " >> ${GNUPLOT_MYDIAG_TMP}
          # All:
          #echo "plot '$MYDIAG_PLOTINFILE_DIFF' using 1:3 with line lw 2 lt rgb \"dark-green\" title \"Diff: ${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\", STATSD_mean lw 2 lc rgb \"gray\" title gprintf( \"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} mean = %g ${MYDIAG_PLOTUDM} \", STATSD_mean ), STATSD_lo_quartile lw 2 lc rgb \"green\" title \"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} 1st quartile\",STATSD_up_quartile lw 2 lc rgb \"green\" title \"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} 3rd quartile\"  " >> ${GNUPLOT_MYDIAG_TMP}
        # SSH
        echo "set ylabel \" Differences ${MYDIAG_PLOTUDM}\"" >> ${GNUPLOT_MYDIAG_TMP}
        echo "set yrange [\"-0.008\":\"0.008\"]" >> ${GNUPLOT_MYDIAG_TMP}
        #echo "plot '$MYDIAG_PLOTINFILE_DIFF' using 1:3 with line lw 2 lt rgb \"${MYDIAG_DIFFCOLOR1}\" title \"${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\", STATSD_mean lw 2 lc rgb \"${MYDIAG_AVG_DIFFCOLOR1}\" title gprintf( \"AVG Diff = %g ${MYDIAG_PLOTUDM} \", STATSD_mean )" >> ${GNUPLOT_MYDIAG_TMP} 
        echo "plot '$MYDIAG_PLOTINFILE_DIFF' using 1:3 with line lw 2 lt rgb \"${MYDIAG_DIFFCOLOR1}\" title \"${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\", STATSD_max lw 0 lc rgb \"${MYDIAG_AVG_DIFFCOLOR1}\" title gprintf( \"MAX Diff = %g [cm/s] \", STATSD_max*100 ), STATSD_mean lw 2 lc rgb \"${MYDIAG_AVG_DIFFCOLOR1}\" title gprintf( \"AVG Diff = %g [cm/s] \", STATSD_mean*100 ), STATSD_min lw 0 lc rgb \"${MYDIAG_AVG_DIFFCOLOR1}\" title gprintf( \"MIN Diff = %g [cm/s] \", STATSD_min*100 ) " >> ${GNUPLOT_MYDIAG_TMP}

        ## Write statistics
        #gnuplot < $GNUPLOT_MYDIAG_TMPTXT || echo "Prob with stat..why?!"
        # Plot
        gnuplot < $GNUPLOT_MYDIAG_TMP >> stat_allvar.txt  || echo "Prob with this plot..why?!"
        #rm -v $GNUPLOT_MYDIAG_TMP
        #rm -v $MYDIAG_FILE_DIFFOUT
       
        else
          echo "ERROR: Mydiag input files NOT found...Why?!"
          echo $MYDIAG_FILE_1 $MYDIAG_FILE_2
        fi
 
     #fi
    IDX_MYDIAG=$(( $IDX_MYDIAG + 1 ))
    done
    
 elif [[ $MYDIAG_FLAG == 2 ]]; then
    echo "Mydiag analysis with 2 datasets.."

    IDX_MYDIAG=0
    for MYDIAG_FIELD in ${MYDIAG_SHORT_NAMES[@]}; do
      #if [[ $IDX_MYDIAG != 2 ]] && [[ $IDX_MYDIAG != 5 ]] && [[ $IDX_MYDIAG != 8 ]] && [[ $IDX_MYDIAG != 11 ]] && [[ $IDX_MYDIAG != 14 ]] && [[ $IDX_MYDIAG != 25 ]] && [[ $IDX_MYDIAG != 33 ]] && [[ $IDX_MYDIAG != 35 ]]; then
        echo "$IDX_MYDIAG ) I am working on ${MYDIAG_SHORT_NAMES[${IDX_MYDIAG}]}.."

        # Infiles
        MYDIAG_FILE_1="${MYDIAG_PATH_1}/${MYDIAG_FILES[${IDX_MYDIAG}]}.nc"
        MYDIAG_FILE_2="${MYDIAG_PATH_2}/${MYDIAG_FILES[${IDX_MYDIAG}]}.nc"
        MYDIAG_FILE_3="${MYDIAG_PATH_3}/${MYDIAG_FILES[${IDX_MYDIAG}]}.nc"
        #     
        if [[ -f $MYDIAG_FILE_1 ]] && [[ -f $MYDIAG_FILE_2 ]] && [[ -f $MYDIAG_FILE_3 ]] ; then

          # Select dates 
          cdo seldate,${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-${ANA_STARTDATE:6:2}T00:00:00,${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}T23:30:00 $MYDIAG_FILE_1 seldate1.nc
          cdo seldate,${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-${ANA_STARTDATE:6:2}T00:00:00,${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}T23:30:00 $MYDIAG_FILE_2 seldate2.nc
          cdo seldate,${ANA_STARTDATE:0:4}-${ANA_STARTDATE:4:2}-${ANA_STARTDATE:6:2}T00:00:00,${ANA_ENDDATE:0:4}-${ANA_ENDDATE:4:2}-${ANA_ENDDATE:6:2}T23:30:00 $MYDIAG_FILE_3 seldate3.nc

          # Compute the differences 
          MYDIAG_DIFFOUT_1="Diff_1_${MYDIAG_FIELD}.nc"
          #cdo sub ${MYDIAG_FILE_1} ${MYDIAG_FILE_2} $MYDIAG_DIFFOUT_1
          cdo sub seldate1.nc seldate2.nc $MYDIAG_DIFFOUT_1
          MYDIAG_DIFFOUT_2="Diff_2_${MYDIAG_FIELD}.nc"
          #cdo sub ${MYDIAG_FILE_3} ${MYDIAG_FILE_2} $MYDIAG_DIFFOUT_2
          cdo sub seldate3.nc seldate2.nc $MYDIAG_DIFFOUT_2
          MYDIAG_DIFFOUT_M="Diff_M_${MYDIAG_FIELD}.nc"
          cdo sub seldate3.nc seldate1.nc $MYDIAG_DIFFOUT_M

          # Extract TS
          MYDIAG_PLOTINFILE_1="plotinfile_1_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_2="plotinfile_2_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_3="plotinfile_3_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_DIFF_1="plotinfile_diff_1_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_DIFF_2="plotinfile_diff_2_${MYDIAG_FIELD}.txt"
          MYDIAG_PLOTINFILE_DIFF_M="plotinfile_diff_M_${MYDIAG_FIELD}.txt"

          # The following two line should be RM 
          #cdo outputtab,date,time,value,name $MYDIAG_FILE_1 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_1
          #cdo outputtab,date,time,value,name $MYDIAG_FILE_2 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_2
          #cdo outputtab,date,time,value,name $MYDIAG_FILE_3 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_3

          cdo outputtab,date,time,value,name seldate1.nc | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_1
          cdo outputtab,date,time,value,name seldate2.nc | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_2
          cdo outputtab,date,time,value,name seldate3.nc | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_3

          cdo outputtab,date,time,value,name $MYDIAG_DIFFOUT_1 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_DIFF_1
          cdo outputtab,date,time,value,name $MYDIAG_DIFFOUT_2 | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_DIFF_2
          cdo outputtab,date,time,value,name $MYDIAG_DIFFOUT_M | grep "$MYDIAG_FIELD" | grep -v "1e+20" > $MYDIAG_PLOTINFILE_DIFF_M
                    

          #ls seldate?.nc
          rm -v seldate1.nc
          rm -v seldate2.nc
          rm -v seldate3.nc


          # Plot TS and Diff
          PLOT_SDATE=$ANA_STARTDATE
          #PLOT_SDATE=$( date -d "${ANA_STARTDATE:0:8} 2 day" +%Y%m%d )
          PLOT_EDATE=$ANA_ENDDATE

          MYDIAG_PLOTUDM=${MYDIAG_UDM[${IDX_MYDIAG}]}

          MYDIAG_PLOT=$( echo "$MYDIAG_PLOT_TPL" | sed -e "s/%FIELD%/${MYDIAG_FILES[${IDX_MYDIAG}]}/g"  -e "s/%DATES%/"${PLOT_SDATE}_${PLOT_EDATE}"/g" )
          echo "MYDIAG_PLOT_TPL $MYDIAG_PLOT"

          GNUPLOT_MYDIAG_TMP="mydiag_tmp.gpl_${MYDIAG_FILES[${IDX_MYDIAG}]}"
          echo "#" > ${GNUPLOT_MYDIAG_TMP}
          #echo "set term jpeg size 1600,1000 giant" >> ${GNUPLOT_MYDIAG_TMP}
          echo "set term jpeg size 1000,800 giant" >> ${GNUPLOT_MYDIAG_TMP}
          echo "set output \"$MYDIAG_PLOT\" " >> ${GNUPLOT_MYDIAG_TMP}

          echo "stats '$MYDIAG_PLOTINFILE_1' using 3 name 'STATS1' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_2' using 3 name 'STATS2' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_3' using 3 name 'STATS3' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_DIFF_1' using 3 name 'STATSD1' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_DIFF_2' using 3 name 'STATSD2' nooutput" >> ${GNUPLOT_MYDIAG_TMP}
          echo "stats '$MYDIAG_PLOTINFILE_DIFF_M' using 3 name 'STATSDM' nooutput" >> ${GNUPLOT_MYDIAG_TMP}

          echo "set multiplot layout 3,1 title \"Time Series ( VAR: ${MYDIAG_LONG_NAMES[${IDX_MYDIAG}]}  DT: $ANA_STARTDATE - $ANA_ENDDATE )\" " >> ${GNUPLOT_MYDIAG_TMP}
          #echo "set title \"Time Series ( VAR: ${MYDIAG_LONG_NAMES[${IDX_MYDIAG}]}  DT: $ANA_STARTDATE - $ANA_ENDDATE )\" " >> ${GNUPLOT_MYDIAG_TMP}
          
          echo "set xlabel \"Time\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xdata time " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set timefmt \"%Y-%m-%d %H:%M:%S\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xrange [\"${PLOT_SDATE:0:4}-${PLOT_SDATE:4:2}-${PLOT_SDATE:6:2} 00:00:00\":\"${PLOT_EDATE:0:4}-${PLOT_EDATE:4:2}-${PLOT_EDATE:6:2} 23:30:00\"] " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set format x \"%m/%Y\" " >> ${GNUPLOT_MYDIAG_TMP} # %d/%m/%Y
          echo "set ylabel \"${MYDIAG_FIELD} ${MYDIAG_PLOTUDM}\" " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set grid " >> ${GNUPLOT_MYDIAG_TMP}
          echo "set xzeroaxis lw 2" >> ${GNUPLOT_MYDIAG_TMP}
          echo "set key right bottom" >> ${GNUPLOT_MYDIAG_TMP} # set key outside

          echo "plot '$MYDIAG_PLOTINFILE_3' using 1:3 with line lw 2 lt rgb \"red\" title gprintf(\"${MYDIAG_INTAG_3} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS3_mean), '$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 2 lt rgb \"orange\" title gprintf(\"${MYDIAG_INTAG_1} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS1_mean), '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 2 lt rgb \"blue\" title gprintf(\"${MYDIAG_INTAG_2} (AVG = %g ${MYDIAG_PLOTUDM} )\", STATS2_mean) " >> ${GNUPLOT_MYDIAG_TMP}
          # SSH:
          #echo "plot '$MYDIAG_PLOTINFILE_3' using 1:3 with line lw 2 lt rgb \"dark-green\" title gprintf(\"${MYDIAG_INTAG_3} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS3_mean), '$MYDIAG_PLOTINFILE_1' using 1:3 with line lw 2 lt rgb \"red\" title gprintf(\"${MYDIAG_INTAG_1} (AVG = %g ${MYDIAG_PLOTUDM} ) \", STATS1_mean), '$MYDIAG_PLOTINFILE_2' using 1:3 with line lw 2 lt rgb \"black\" title gprintf(\"${MYDIAG_INTAG_2} (AVG = %g ${MYDIAG_PLOTUDM} )\", STATS2_mean) " >> ${GNUPLOT_MYDIAG_TMP}
          ##echo "plot '$MYDIAG_PLOTINFILE_DIFF_1' using 1:3 with line lw 2 lt rgb \"black\" title \"Diff: ${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\",'$MYDIAG_PLOTINFILE_DIFF_2' using 1:3 with line lw 2 lt rgb \"green\" title \"Diff: ${MYDIAG_INTAG_3} - ${MYDIAG_INTAG_2}\" " >> ${GNUPLOT_MYDIAG_TMP}
          #echo "set yrange [\"-0.10\":\"0.10\"] " >> ${GNUPLOT_MYDIAG_TMP}
          echo "plot '$MYDIAG_PLOTINFILE_DIFF_1' using 1:3 with line lw 2 lt rgb \"orange\" title \"Diff: ${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\",'$MYDIAG_PLOTINFILE_DIFF_2' using 1:3 with line lw 2 lt rgb \"red\" title \"Diff: ${MYDIAG_INTAG_3} - ${MYDIAG_INTAG_2}\"" >> ${GNUPLOT_MYDIAG_TMP}
          echo "plot '$MYDIAG_PLOTINFILE_DIFF_M' using 1:3 with line lw 2 lt rgb \"dark-green\" title \"Diff: ${MYDIAG_INTAG_3} - ${MYDIAG_INTAG_1}\"" >> ${GNUPLOT_MYDIAG_TMP}
          # OLD with mean diff and quartiles:
          #echo "plot '$MYDIAG_PLOTINFILE_DIFF_1' using 1:3 with line lw 2 lt rgb \"orange\" title \"Diff: ${MYDIAG_INTAG_1} - ${MYDIAG_INTAG_2}\",'$MYDIAG_PLOTINFILE_DIFF_2' using 1:3 with line lw 2 lt rgb \"red\" title \"Diff: ${MYDIAG_INTAG_3} - ${MYDIAG_INTAG_2}\", STATSD1_mean lw 2 lc rgb \"black\" title gprintf(\"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} mean %g ${MYDIAG_PLOTUDM} \", STATSD1_mean) , STATSD1_lo_quartile lw 2 lc rgb \"orange\" title \"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} 1st quartile\",STATSD1_up_quartile lw 2 lc rgb \"orange\" title \"${MYDIAG_INTAG_1}-${MYDIAG_INTAG_2} 3rd quartile\",STATSD2_mean lw 2 lc rgb \"blue\" title gprintf( \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_2} mean = %g ${MYDIAG_PLOTUDM} \", STATSD2_mean ), STATSD2_lo_quartile lw 2 lc rgb \"red\" title \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_2} 1st quartile\",STATSD2_up_quartile lw 2 lc rgb \"red\" title \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_2} 3rd quartile\"  " >> ${GNUPLOT_MYDIAG_TMP}
         #echo "plot '$MYDIAG_PLOTINFILE_DIFF_M' using 1:3 with points pt 1 ps 1 lt rgb \"dark-green\" title \"Diff: ${MYDIAG_INTAG_3} - ${MYDIAG_INTAG_1}\",STATSDM_mean lw 2 lc rgb \"black\" title gprintf( \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_1} mean = %g ${MYDIAG_PLOTUDM} \", STATSDM_mean ), STATSDM_lo_quartile lw 2 lc rgb \"green\" title \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_1} 1st quartile\",STATSDM_up_quartile lw 2 lc rgb \"green\" title \"${MYDIAG_INTAG_3}-${MYDIAG_INTAG_1} 3rd quartile\" " >> ${GNUPLOT_MYDIAG_TMP}

        # Plot
        gnuplot < $GNUPLOT_MYDIAG_TMP || echo "Prob with this plot..why?!"
        #rm -v $GNUPLOT_MYDIAG_TMP
        #rm -v $MYDIAG_FILE_DIFFOUT

        else
          echo "ERROR: Mydiag input files NOT found...Why?!"
          echo $MYDIAG_FILE_1 $MYDIAG_FILE_2 $MYDIAG_FILE_3
        fi

     #fi
    IDX_MYDIAG=$(( $IDX_MYDIAG + 1 ))
    done


 fi

 echo "WORKDIR: $ANA_WORKDIR"

###################### POSTPROC ###########################

# Output check

# Archive


