#!/bin/bash
#
# ACG 19/01/2020
# Script for detiding with Doodson filter 
# Ini file: filter.ini 
#
#set -u
set -e
#set -x 
################### ENV SETTINGS ##############################
SRC_DIR="/users/home/ag15419/tides_pp/doodson_filter"
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
         rm $ANA_WORKDIR/$TO_BE_RM
         echo $TO_BE_RM
     done
     # Cp the ini file to the workdir 
     cp $INIFILE ${ANA_WORKDIR}

  else
     echo "ERROR: WORKDIR $ANA_WORKDIR NOT FOUND!!"
     exit
  fi

  # Input File check and link

    if [[ -d ${ANA_INPATH} ]]; then

      IDX_DATE=$ANA_STARTDATE
      
      while [[ $IDX_DATE -le $ANA_ENDDATE ]]; do

        echo "Date: $IDX_DATE"
        ANA_INFILE=$( echo ${ANA_INFILE_TPL} | sed -e s/%YYYYMMDD%/$IDX_DATE/g  )

        FOUND_NUM=0
        if [[ $GRID_TO_EXTRACT != "uv2t" ]]; then
           if [[ -e ${ANA_INPATH}/${IDX_DATE:0:6}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATH}/${IDX_DATE:0:6}/$ANA_INFILE .
           else
              echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATH}/${IDX_DATE:0:6}/"
           fi
        else
           if [[ -e ${ANA_INPATH}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATH}/$ANA_INFILE .
           fi
        fi

      IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d ) 
      done    
 
    else 
      echo "ERROR: Input dir ${ANA_INPATH}]} NOT FOUND!!"
      exit
    fi


    # Set the environment
  
    module load $DOODF_MODULE
    echo "Loading the Enviroment.."

#############################################################################
####################### DOODSON FILTER ######################################

 echo "------ Applying Doodson filter to model outputs ------"
 echo "The Filter has a window of ${#F_NUM[@]} hours before and after 12:00 "
 echo "The values of the filter are [Zanella, 2015]: ${F_NUM[@]}"

############################################################################

 if [[ ${#F_NUM[@]} -lt 24 ]]; then

    echo "I am working on ${ANA_INTAG} dataset..."

    echo "I am working on $VAR_NAME field..."

    # Input file names (linked in work dir)
    EXT_INS=$(  echo "${ANA_INFILE_TPL}" | sed -e "s/%YYYYMMDD%/"*"/g" )

    VLEV=$VAR_LEV
    echo "Depth: $VLEV .. "

    echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"
    
    # Field selection and Hours splitting 
 
      # Loop on input daily files
      IDX_NC=0
      for NF_SPLIT in $EXT_INS; do
          echo "Infile: $NF_SPLIT"

          # Select var and split hours
          echo "Field selection and hours splitting..in ${NF_SPLIT}_ files.."
          cdo selname,${VAR_NAME} $NF_SPLIT ${NF_SPLIT}_${VAR_NAME} 
 
          cdo splithour ${NF_SPLIT}_${VAR_NAME} ${NF_SPLIT}_
          rm -v ${NF_SPLIT}_${VAR_NAME}

      IDX_NC=$(( $IDX_NC + 1 ))
      done
 
    # 0) Loop on input files
    # IDX initialization (since the Doodson filter covers more than )
    IDX_DATE=$( date -d "$ANA_STARTDATE 1 day" +%Y%m%d )

    while [[ $IDX_DATE -lt $ANA_ENDDATE ]]; do

          # Doodson Filter outfile name 
          FILTER_OUTFILE=$( echo "$FILTER_OUTFILE_TPL" | sed -e "s/%FIELD%/${VAR_NAME}/g" -e "s/%INDATASET%/${ANA_INTAG}/g"  -e "s/%YYYYMMDD%/${IDX_DATE}/g"  )
          echo "I am working on day: ${IDX_DATE}"
          echo "The outfile storing the Doddson filtered values is: ${FILTER_OUTFILE}"

          # 1) Loop on hours (involved in the filter computation)
          H_IDX=0

          # Element of the doodson filter sum   
          SUM_EL=1 

          # Num of hours outside the current day
          MAX_HOUT=$(( ${#F_NUM[@]} - 12 )) 
          MIN_HOUT=$(( 24 - $MAX_HOUT ))
          PREV_HOUR=$MIN_HOUT
          SUCC_HOUR=$MAX_HOUT
 

          # Hours outside the current day
          echo "###################"  
          while [[ $H_IDX -lt $MAX_HOUT  ]]; do
                echo "===== ====="
                echo "Computing the ${H_IDX} element of Doodson filter (outside the current day).."
                
                PREV_DAY=$( date -d "$IDX_DATE -1 day" +%Y%m%d )
                PREV_HOUR=$( printf %02d $PREV_HOUR )
                FILE_LOWVAL=$(  echo "${ANA_INFILE_TPL}_${PREV_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$PREV_DAY/g" )

                SUCC_DAY=$( date -d "$IDX_DATE 1 day" +%Y%m%d )
                SUCC_HOUR=$( printf %02d $SUCC_HOUR )
                FILE_HIVAL=$(  echo "${ANA_INFILE_TPL}_${SUCC_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$SUCC_DAY/g" )
                
                # Multiply each extracted hour for F doodson num
                echo "Working on $PREV_DAY ${PREV_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_LOWVAL tmp_Lel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))
                echo "Working on $SUCC_DAY ${SUCC_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_HIVAL tmp_Hel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))

          PREV_HOUR=$(( $PREV_HOUR + 1 ))
          SUCC_HOUR=$(( $SUCC_HOUR - 1 ))
          H_IDX=$(( $H_IDX + 1 ))
          done

          # 00 Case (Needs a specific treatment because the 24 are the 00 of the succ day!)
          echo "###################"
          if [[ $H_IDX == 7 ]]; then
             echo "===== ====="
             echo "Computing the ${H_IDX} element of Doodson filter (00 case).."
                PREV_DAY=$IDX_DATE
                PREV_HOUR="00"
                FILE_LOWVAL=$(  echo "${ANA_INFILE_TPL}_${PREV_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$PREV_DAY/g" )

                SUCC_DAY=$( date -d "$IDX_DATE 1 day" +%Y%m%d )
                SUCC_HOUR="00"
                FILE_HIVAL=$(  echo "${ANA_INFILE_TPL}_${SUCC_HOUR}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/$SUCC_DAY/g" )

                # Multiply each extracted hour for F doodson num
                echo "Working on $PREV_DAY ${PREV_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_LOWVAL tmp_Lel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))
                echo "Working on $SUCC_DAY ${SUCC_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_HIVAL tmp_Hel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))

          H_IDX=$(( $H_IDX + 1 ))
          fi
          # Hours inside the current day          
          echo "###################"
          PREV_HOUR=1
          SUCC_HOUR=23
          while [[ $H_IDX -lt ${#F_NUM[@]}  ]]; do
                echo "===== ====="
                echo "Computing the ${H_IDX} element of Doodson filter (inside the current day).."
          
                PREV_HOUR_2=$( printf %02d $PREV_HOUR )
                FILE_LOWVAL=$(  echo "${ANA_INFILE_TPL}_${PREV_HOUR_2}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/${IDX_DATE}/g" )

                SUCC_HOUR_2=$( printf %02d $SUCC_HOUR )
                FILE_HIVAL=$(  echo "${ANA_INFILE_TPL}_${SUCC_HOUR_2}.nc${VAR_NAME}" | sed -e "s/%YYYYMMDD%/${IDX_DATE}/g" )

                # Multiply each extracted hour for F doodson num
                echo "Working on ${IDX_DATE} ${PREV_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_LOWVAL tmp_Lel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))
                echo "Working on ${IDX_DATE} ${SUCC_HOUR} (Doodson num ${F_NUM[$H_IDX]})"
                cdo -mulc,${F_NUM[$H_IDX]} $FILE_HIVAL tmp_Hel_${SUM_EL}.nc
                SUM_EL=$(( $SUM_EL + 1 ))

          PREV_HOUR=$(( $PREV_HOUR + 1 ))
          SUCC_HOUR=$(( $SUCC_HOUR - 1 ))
          H_IDX=$(( $H_IDX + 1 ))
          done

          echo "Computing the sum of elements.."
          # 2) Cat all elements
          cdo mergetime tmp_*.nc tot_el.nc

          # 3) sum elements and multiply for 1/30 (as prescrived by Doodson filter formulation)
          cdo timsum tot_el.nc sumtot.nc
          cdo divc,30.0 sumtot.nc ${FILTER_OUTFILE}
          #cdo -divc,30.0 -timsum tot_el.nc ${FILTER_OUTFILE}

          # Clean the workdir (must be done before the next daily cicle..)
          rm -v tmp_*.nc
          rm -v tot_el.nc
          rm -v sumtot.nc

     IDX_DATE=$( date -d "$IDX_DATE 1 day" +%Y%m%d )
     done

 fi




###################### POSTPROC ###########################

# Output check

# Archive


