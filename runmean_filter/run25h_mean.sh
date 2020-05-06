#!/bin/bash
#
# ACG 19/01/2020
# Script for detiding with 25h running mean
# Ini file: filter.ini 
#
#set -u
set -e
#set -x 
################### ENV SETTINGS ##############################
SRC_DIR="/users/home/ag15419/tides_pp/runmean_filter"
################### PREPROC ###################################

# Source ini file
  INIFILE="${SRC_DIR}/filter.ini"
  source $INIFILE
  echo "source $INIFILE ... Done!"

# Read and check infos (work dir, file names, archive dir, etc.)

  # Workdir check
  if [[ -d $ANA_WORKDIR ]]; then
     cd $ANA_WORKDIR
     echo "WORKDIR: $ANA_WORKDIR"
     
     # Clean workdir
     echo "WARNING: I am going to remove all files in $ANA_WORKDIR ..."
     sleep 10
     for TO_BE_RM in $( ls $ANA_WORKDIR ); do
         #rm $ANA_WORKDIR/$TO_BE_RM
         echo $TO_BE_RM
     done
     # Cp the ini file to the workdir 
     cp $INIFILE ${ANA_WORKDIR}

  else
     echo "ERROR: WORKDIR $ANA_WORKDIR NOT FOUND!!"
     exit
  fi

  # Input File check and link

#    if [[ -d ${ANA_INPATH} ]]; then
#
#      IDX_DATE=$ANA_STARTDATE
#      
#      while [[ $IDX_DATE -le $ANA_ENDDATE ]]; do
#
#        echo "Date: $IDX_DATE"
#        ANA_INFILE=$( echo ${ANA_INFILE_TPL} | sed -e s/%YYYYMMDD%/$IDX_DATE/g  )
#
#        FOUND_NUM=0
#        if [[ $GRID_TO_EXTRACT != "uv2t" ]]; then
#           if [[ -e ${ANA_INPATH}/${IDX_DATE:0:6}/$ANA_INFILE ]]; then
#           #if [[ -e ${ANA_INPATH}/$ANA_INFILE ]]; then
#              FOUND_NUM=$(( $FOUND_NUM + 1 ))
#              echo "Found infile: $ANA_INFILE"
#              ln -sf ${ANA_INPATH}/${IDX_DATE:0:6}/$ANA_INFILE .
#           else
#              echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATH}/${IDX_DATE:0:6}/"
#           fi
#        else
#           if [[ -e ${ANA_INPATH}/$ANA_INFILE ]]; then
#              FOUND_NUM=$(( $FOUND_NUM + 1 ))
#              echo "Found infile: $ANA_INFILE"
#              ln -sf ${ANA_INPATH}/$ANA_INFILE .
#           fi
#        fi
#
#      IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d ) 
#      done    
# 
#    else 
#      echo "ERROR: Input dir ${ANA_INPATH}]} NOT FOUND!!"
#      exit
#    fi


    # Set the environment
  
    module load $DOODF_MODULE
    echo "Loading the Enviroment.."

#############################################################################
####################### FILTER ######################################

 echo "------ Applying HOURLY 25h (up to 24 or lower!) Running-mean filter to model outputs ------"
 echo "The Filter has a window of ${#HF_NUM[@]} hours before and after each hour "
 echo "The values of the filter are: ${HF_NUM[@]}"

############################################################################

 if [[ ${#HF_NUM[@]} -lt 25 ]]; then

    echo "I am working on ${ANA_INTAG} dataset..."

    echo "I am working on $VAR_NAME field..."

    # Input file names (linked in work dir)
    EXT_INS=$(  echo "${ANA_INFILE_TPL}" | sed -e "s/%YYYYMMDD%/"*"/g" )

    VLEV=$VAR_LEV
    echo "Depth: $VLEV .. "

    echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"
    
#    # Field selection and Hours splitting 
# 
#      # Loop on input daily files
#      IDX_NC=0
#      for NF_SPLIT in $EXT_INS; do
#          echo "Infile: $NF_SPLIT"
#
#          # Select var and split hours
#          echo "Field selection and hours splitting..in ${NF_SPLIT}_ files.."
#          cdo selname,${VAR_NAME} $NF_SPLIT ${NF_SPLIT}_${VAR_NAME} 
# 
#          cdo splithour ${NF_SPLIT}_${VAR_NAME} ${NF_SPLIT}_
#          rm -v ${NF_SPLIT}_${VAR_NAME}
#
#      IDX_NC=$(( $IDX_NC + 1 ))
#      done
 
    # Loop on input files
    # IDX initialization (since the 25h running mean filter covers 5 days)
    IDX_DATE=$( date -d "$ANA_STARTDATE 1 day" +%Y%m%d )
    IDXEND_DATE=$( date -d "$ANA_ENDDATE -1 day" +%Y%m%d )

    while [[ $IDX_DATE -le $IDXEND_DATE ]]; do


      # Running-mean Filter outfile name 
      HFILTER_OUTFILE=$( echo "$HFILTER_OUTFILE_TPL" | sed -e "s/%FIELD%/${VAR_NAME}/g" -e "s/%INDATASET%/${ANA_INTAG}/g"  -e "s/%YYYYMMDD%/${IDX_DATE}/g" )
      echo "I am working on day: ${IDX_DATE}"
      echo "The outfile storing the Running-mean filtered values is: ${HFILTER_OUTFILE}"

      # Loop on HOURS 
      DF_HOUR=0
      DF_LAST_HOUR=24
      while [[ $DF_HOUR -lt $DF_LAST_HOUR ]]; do

            DF_HOUR_2DIG=$( printf %02d $DF_HOUR )
            echo "I am working on hour: ${DF_HOUR_2DIG}"

            # Loop on filter elements 
            H_IDX=00
            while [[ $H_IDX -lt ${#HF_NUM[@]} ]]; do
           
                DELTA_H=$(( $H_IDX + 1 ))                

                PREV_DATE=$( date -d "${IDX_DATE} ${DF_HOUR_2DIG} -${DELTA_H} hour" +%Y%m%d%H )
                PREV_DAY=${PREV_DATE:0:8}
                PREV_HOUR=${PREV_DATE:8:2}
                FILE_PREVVAL=$(  echo "${ANA_INFILE_TPL}_${PREV_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$PREV_DAY/g" )

                SUCC_DATE=$( date -d "${IDX_DATE} ${DF_HOUR_2DIG} ${DELTA_H} hour" +%Y%m%d%H )
                SUCC_DAY=${SUCC_DATE:0:8}
                SUCC_HOUR=${SUCC_DATE:8:2}
                FILE_SUCCVAL=$(  echo "${ANA_INFILE_TPL}_${SUCC_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$SUCC_DAY/g" )


                # Multiply each extracted hour for F doodson num
                echo "Working on $PREV_DAY ${PREV_HOUR} (Coef. ${HF_NUM[$H_IDX]})"
                cdo -mulc,${HF_NUM[$H_IDX]} $FILE_PREVVAL tmp_${H_IDX}_Pre.nc
                echo "Working on $SUCC_DAY ${SUCC_HOUR} (Coef. ${HF_NUM[$H_IDX]})"
                cdo -mulc,${HF_NUM[$H_IDX]} $FILE_SUCCVAL tmp_${H_IDX}_Suc.nc


            H_IDX=$(( $H_IDX + 1 ))
            done

            echo "Computing the sum of elements (day: ${IDX_DATE} hour: ${DF_HOUR_2DIG} ) .."
            # Cat all elements
            cdo mergetime tmp_*.nc tot_el.nc

            # Sum elements and multiply for 1/(2*h+1) (as prescrived by Running-mean filters formulation)
            DIV_NUM=$(( ${#HF_NUM[@]} * 2 + 1 ))
            echo "(2*h+1) = $DIV_NUM"
            cdo timsum tot_el.nc sumtot.nc
            cdo divc,$DIV_NUM sumtot.nc det_${IDX_DATE}_${DF_HOUR_2DIG}.nc

            # Clean the workdir (must be done before the next daily cicle..)
            rm -v tmp_*.nc
            rm -v tot_el.nc
            rm -v sumtot.nc

      DF_HOUR=$(( $DF_HOUR + 1 ))
      done
      
      # Cat hourly values in a single daily file 
      cdo mergetime det_${IDX_DATE}_*.nc ${HFILTER_OUTFILE}

      # Clean the workdir (must be done before the next daily cicle..)
      rm -v det_${IDX_DATE}_*.nc

 
     IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d )
     done

 fi




###################### POSTPROC ###########################

# Output check

# Archive


