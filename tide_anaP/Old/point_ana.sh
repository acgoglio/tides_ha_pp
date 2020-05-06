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

  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     
     # Clean workdir
     echo "WARNING: I am going to remove all files in $ANA_WORKDIR ..."
     sleep 10
     for TO_BE_RM in $( ls $ANA_WORKDIR ); do
         rm -v $ANA_WORKDIR/$TO_BE_RM
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
           FOUND_NUM=$(( $FOUND_NUM + 1 ))
           echo "Found infile: $ANA_INFILE"
           ln -sf ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE .
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
  
  if [[ $DMMM_FLAG == 1 ]]; then
    echo "Min, mean and max analisys over the whole domain is required!"
    module load $DMMM_MODULE
  fi

  if [[ $TS_FLAG == 1 ]]; then
    echo "TS analisys for points listed in $TS_COOFILE is required!"
    module load $TS_MODULE
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
      echo "Prova EXT_INS: $EXT_INS"
      DMMM_OUTFILE=$( echo "$DMMM_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g"  )
      echo "Prova $DMMM_OUTFILE"

      echo "# DMMM extraction form ${ANA_INTAG[$IDX_IN]}" > $DMMM_OUTFILE 
      echo "# " >> $DMMM_OUTFILE
      echo "# -1 :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name" >> $DMMM_OUTFILE

      # Extraction
      for DMMEXT in $EXT_INS; do   
        echo "cdo infon $DMMEXT | grep -v "Date" | grep "$VAR"  >> $DMMM_OUTFILE"
        cdo infon $DMMEXT | grep -v "Date" | grep "$VAR"  >> $DMMM_OUTFILE
        rm -v $DMMEXT
      done

     done

    IDX_IN=$(( $IDX_IN + 1 ))
    done 

    #########################
    echo "------ Plot ------"
    echo "PROV ${DMMM_FIELDS[@]}"
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

        GNUPLOT_TMP="tmp.gpl"

        echo "PROVA: 01 $TXT_TAB_0 02 $TXT_TAB_02 $DMMM_PLOT 2 $PVAR 3 $ANA_STARTDATE 4 $ANA_ENDDATE 5 $PLOT_SDATE 6 $PLOT_EDATE 7 $PLOT_VAR 8 ${PLOT_UNITS}"

#        cat << EOF > ${GNUPLOT_TMP}  

        echo "#" > ${GNUPLOT_TMP}
        echo "set term jpeg size 1000,600" >> ${GNUPLOT_TMP}
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
          echo "Prova: sto definendo i valori di indice ${COO_IDX} "
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

    echo "------ TS Extraction ------"

  TS_OUTFILE_PRE="ts"
  TS_OUTFILE_TPL="${TS_OUTFILE_PRE}_${ANA_STARTDATE}_${ANA_ENDDATE}_%STZ%_%FIELD%_%IDX%.txt"

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    echo "I am workin on ${ANA_INTAG[$IDX_IN]} dataset..."

    for VAR in ${TS_FIELDS[@]}; do

     echo "Extracting field: $VAR .. "

     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )
     echo "Prova $EXT_INS"
     
     TS_IDX=0
     while [[ $TS_IDX -lt $COO_IDX ]]; do

      TS_OUTFILE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/$IDX_IN/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
      echo "Prova $TS_OUTFILE"

      echo "# TS extraction from ${ANA_INTAG[$IDX_IN]}" > $TS_OUTFILE
      echo "# " >> $TS_OUTFILE
      echo "# STZ: ${COO_TS_NAME[$TS_IDX]}" >> $TS_OUTFILE
      echo "# VAR: ${VAR}" >> $TS_OUTFILE
      echo "# " >> $TS_OUTFILE
      echo "# -1 :       Date     Time   Level Gridsize    Miss :     Minimum        Mean     Maximum : Parameter name" >> $TS_OUTFILE

      # Extraction
      for TSEXT in $EXT_INS; do
        
        # Time interpolation (for comparison with obs)
        cdo intntime,2 $TSEXT 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc 

        # With land/sea mask (Nan)
        cdo setctomiss,0.0000 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        ###cdo setctomiss,0.0000 $TSEXT miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo -remapdis,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        # Mean comp
        cdo timmean tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        cdo sub tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        #rm -v 00_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        #rm -v miss_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        #rm -v mean_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc
        #rm -v tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        ####cdo -infon zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc  | grep "$VAR" | grep "00:00 " >> $TS_OUTFILE
        cdo -infon zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc  | grep "$VAR" >> $TS_OUTFILE 
        #rm -v zero_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        # Without mask:
        #echo "PROVA: cdo -remapnn,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} $TSEXT tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc" 
        #cdo -remapnn,lon=${COO_TS_LON[$TS_IDX]}/lat=${COO_TS_LAT[$TS_IDX]} $TSEXT tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc

        #cdo -infon tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc  | grep "$VAR" >> $TS_OUTFILE
        #rm -v tsext_tmp_${TS_IDX}_${COO_TS_NAME[$TS_IDX]}.nc


      done

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
	     echo "Prova: cdo infon $OBS_FILE | grep "$O_VAR" >> $OBS_OUTFILE "


             # with zero sub
             for OBS_FOUND in $( ls $OBS_FILE ); do
              if [[ -f $OBS_FOUND ]]; then
                cdo timmean $OBS_FILE obs_mean.nc
                cdo sub $OBS_FILE obs_mean.nc obs_zero.nc
                cdo infon obs_zero.nc | grep "${O_VAR} " >> $OBS_OUTFILE
                rm -v obs_mean.nc
                rm -v obs_zero.nc
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
             echo "Prova: OBS FILE = $OBS_FILE "

             for OBS_FOUND in $( ls $OBS_FILE ); do
              if [[ -f $OBS_FOUND ]]; then
                while read ISPRA_LINE ; do
                   if [[ ${ISPRA_LINE:0:1} != "G" ]]; then
                      ISPRA_DATA=${ISPRA_LINE:0:8} 
                      ISPRA_ORA=$( echo $ISPRA_LINE | cut -f 2 -d";" )
                      ISPRA_VAL=$( echo $ISPRA_LINE | cut -f 3 -d";" )
                      #ISPRA_VAL_M=$( echo "${ISPRA_VAL:0:4} / 10 " | bc -l )
                      echo " $ID_LINE : 20${ISPRA_DATA:6:2}-${ISPRA_DATA:3:2}-${ISPRA_DATA:0:2} ${ISPRA_ORA}:00       0        1       0 :                 ${ISPRA_VAL}    : sossheig " >> $OBS_OUTFILE
                   
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

# Plot ts

  echo "------ Plot TS ------"

    echo "I am workin on # ${INSET_NUM} dataset(s)..."

    ID_V=0
    for VAR in ${TS_FIELDS[@]}; do
      O_VAR=${OBS_VAR[${ID_V}]}
      ID_V=$(( $ID_V + 1 ))  
      echo "Prova coppia var: $VAR e $O_VAR"

      IDX_UDMVAR=0
      echo "Plot field: $VAR .. "

      PLOT_UNITS=${DMMM_UDM[$IDX_PVAR]}
      PLOT_SDATE=$ANA_STARTDATE
      PLOT_EDATE=$ANA_ENDDATE


      TS_IDX=0
      echo "Prova: COO_IDX=$COO_IDX"
      while [[ $TS_IDX -lt $COO_IDX ]]; do
      echo "STZ: ${COO_TS_NAME[$TS_IDX]}"

        TSPLOT_INFILE_0=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/0/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
        TSPLOT_INFILE_1=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/1/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
        TSPLOT_INFILE_2=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/2/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
 
        if [[ $DIFF_FLAG == 1 ]] && [[ $OBS_FLAG == 0 ]]; then
          echo "DIFF active"
          INFILE_PASTE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/all/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
            paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_1 > $INFILE_PASTE
            echo "INFILE_PASTE=$INFILE_PASTE"         
        fi

        if [[ $OBS_FLAG != 0 ]]; then
           TSPLOT_INFILE_OBS=$( echo "$OBS_OUTFILE_TPL" | sed -e "s/%FIELD%/${O_VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
           if [[ $DIFF_FLAG == 1 ]]; then
              INFILE_PASTE=$( echo "$TS_OUTFILE_TPL" | sed -e "s/%IDX%/allobs/g" -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" )
              if [[ ${MOD1_FLAG} == 1 ]]; then
                 paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_OBS > $INFILE_PASTE
              else
                 paste $TSPLOT_INFILE_0 $TSPLOT_INFILE_1 $TSPLOT_INFILE_OBS > $INFILE_PASTE
              fi
           fi
        fi

        echo "Prova $TSPLOT_INFILE_0"
        echo "Prova $TSPLOT_INFILE_1"

        TS_PLOT=$( echo "$TS_PLOT_TPL" | sed -e "s/%FIELD%/${VAR}/g" -e "s/%STZ%/${COO_TS_NAME[$TS_IDX]}/g" -e "s/%DATES%/"${PLOT_SDATE}_${PLOT_EDATE}"/g" )
       
        if [[ $DIFF_FLAG == 0 ]]; then 

          # Plot mod ts 
          GNUPLOT_TS_TMP="ts_tmp.gpl"    

          echo "#" > ${GNUPLOT_TS_TMP}
          echo "set term jpeg size 1500,600" >> ${GNUPLOT_TS_TMP}
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
             #3 mod
             if [[ ${#ANA_INTAG[@]} == 3 ]]; then
               echo "3 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_2' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[2]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100 - STATS_mean/100 )  with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}

             #2 mod
             elif [[ ${#ANA_INTAG[@]} == 2 ]]; then
               echo "2 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"orange\" title \"${ANA_INTAG[1]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100 - STATS_mean/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}

             # 1 model..
             elif [[ ${#ANA_INTAG[@]} == 1 ]]; then
               echo "1 mod"
               echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"red\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_OBS' using 3:(\$9/100-STATS_mean/100) with line lw 2 lt rgb \"blue\" title \"Obs\" " >> ${GNUPLOT_TS_TMP}
             fi
           
          fi

        else 
          # + Diff

          # Plot mod ts 
          GNUPLOT_TS_TMP="ts_tmp.gpl"

          echo "#" > ${GNUPLOT_TS_TMP}
          echo "set term jpeg size 1500,900" >> ${GNUPLOT_TS_TMP}
          echo "set output \"$TS_PLOT\" " >> ${GNUPLOT_TS_TMP}
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
             echo "plot '$TSPLOT_INFILE_0' using 3:9 with line lw 2 title \"Diff: ${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
             echo "plot '$INFILE_PASTE' using 3:(\$9-\$20) with line lw 2 title \"${ANA_INTAG[0]}-${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
          elif [[ $OBS_FLAG != 0 ]]; then
             echo "Obs plot"
             echo "plot '$TSPLOT_INFILE_OBS' using 3:(0):9 with filledc lt rgb \"gray50\" title \"Obs\", '$TSPLOT_INFILE_0' using 3:9 with line lw 2 lt rgb \"black\" title \"${ANA_INTAG[0]}\", '$TSPLOT_INFILE_1' using 3:9 with line lw 2 lt rgb \"green\" title \"${ANA_INTAG[1]}\" " >> ${GNUPLOT_TS_TMP}
             #echo "plot '$INFILE_PASTE' using 3:(\$9-\$31) with line lw 2 title \"${ANA_INTAG[0]}-OBS\", '$INFILE_PASTE' using 3:(\$20-\$31) with line lw 2 title \"${ANA_INTAG[1]}-OBS\" " >> ${GNUPLOT_TS_TMP}
             #echo "plot '$INFILE_PASTE' using 3:(0):(\$9-\$31) with filledc lc rgb \"red\" title \"${ANA_INTAG[0]}-OBS\", '$INFILE_PASTE' using 3:(0):(\$20-\$31) with filledc lc rgb \"blue\" title \"${ANA_INTAG[1]}-OBS\" " >> ${GNUPLOT_TS_TMP}
             #echo "plot '$INFILE_PASTE' using 3:(\$9-\$31):(\$20-\$31) with filledc lc rgb \"gray90\" title \"${ANA_INTAG[0]}-${ANA_INTAG[1]}\", '$INFILE_PASTE' using 3:(\$9-\$31) with line lw 2 lc rgb \"red\" title \"${ANA_INTAG[0]}-OBS\",'$INFILE_PASTE' using 3:(\$20-\$31) with line lw 2 lc rgb \"blue\" title \"${ANA_INTAG[1]}-OBS\" " >> ${GNUPLOT_TS_TMP}
             echo "plot '$INFILE_PASTE' using 3:(\$9-\$31) with line lw 2 lc rgb \"red\" title \"${ANA_INTAG[0]}-OBS\",'$INFILE_PASTE' using 3:(\$20-\$31) with line lw 2 lc rgb \"orange\" title \"${ANA_INTAG[1]}-OBS\" " >> ${GNUPLOT_TS_TMP}
             #echo "plot '$TSPLOT_INFILE_OBS' using 3:9 with line lw 2 title \"Obs\" " >> ${GNUPLOT_TS_TMP}
          fi
 
        fi

        # Plot
        gnuplot < $GNUPLOT_TS_TMP || echo "Prob with this plot..why?!"
        #rm -v $GNUPLOT_TS_TMP

      TS_IDX=$(( $TS_IDX + 1 ))
      done

     IDX_UDMVAR=$(( $IDX_UDMVAR + 1 ))
     done

 fi


# +-----------------------+
# | MYDIAG DIFF           |
# +-----------------------+ 

 #if [[]]; then

 #fi


###################### POSTPROC ###########################

# Output check

# Archive


