#!/bin/bash
#
# ACG 06/11/2019
# Script for sections extraction 
# Ini file: sect_extr.ini 
#
#set -u
set -e
#set -x 
################### PREPROC ###################################

# Source ini file
  SRC_DIR="/users/home/ag15419/tides_pp/sections"
  INIFILE="${SRC_DIR}/sect_extr.ini"
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

     # Cp ini file in the workdir
     cp $INIFILE ${ANA_WORKDIR}

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

  # Input File check and link

  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    if [[ -d ${ANA_INPATHS[${IDX_IN}]} ]]; then

      IDX_DATE=$ANA_STARTDATE
      
      while [[ $IDX_DATE -le $ANA_ENDDATE ]]; do

        echo "Date: $IDX_DATE"
        ANA_INFILE=$( echo ${ANA_INFILES_TPL[${IDX_IN}]} | sed -e s/%YYYYMMDD%/$IDX_DATE/g  )

        FOUND_NUM=0
        if [[ $GRID_TO_EXTRACT != 'uv2t' ]]; then
           if [[ -e ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/$ANA_INFILE .
           else
              echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATHS[$IDX_IN]}/${IDX_DATE:0:6}/"
           fi
        else
           if [[ -e ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE ]]; then
              FOUND_NUM=$(( $FOUND_NUM + 1 ))
              echo "Found infile: $ANA_INFILE"
              ln -sf ${ANA_INPATHS[$IDX_IN]}/$ANA_INFILE .
           else
              echo "NOT Found infile: $ANA_INFILE in path: ${ANA_INPATHS[$IDX_IN]}/"
           fi
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
  
  if [[ $SECT_FLAG == 1 ]]; then
    echo "Map analisys required!"
    module load $SECT_MODULE
    echo "Loading the Enviroment.."
  fi


####################### ANALISYS ##############################

# Extract values from netCDF

 if [[ $SECT_FLAG == 1 ]]; then
  
  echo "LAT/LON Reading..."
  echo "LAT RANGE = ${LAT_RANGE[@]}"
  echo "LON RANGE = ${LON_RANGE[@]}"
  

  echo "------ Section Extraction from model outputs ------"

  # Loop on model datasets
  IDX_IN=0
  while [[ $IDX_IN -lt $INSET_NUM ]]; do

    echo "I am working on ${ANA_INTAG[$IDX_IN]} dataset..."

    ############# 3D VARS #########################

    echo "I am working on 3D vars..."

    # Loop on 3D FIELDS
    for VAR in ${VAR3D_NAME[@]}; do
     echo "Extracting field: $VAR .. "
     
     # Input file names (linked in work dir)
     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )
     
     SECT_OUTFILE=$( echo "$SECT_OUTFILE_TPL" | sed  -e "s/%SECT_LABEL%/${SECT_LABEL}/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

     echo "# SECTION extraction from ${ANA_INTAG[$IDX_IN]}" 
     echo "# SECT_LABEL: ${SECT_LABEL}" 
     echo "# VAR: ${VAR}" 
     echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"
     echo "# Seasonal analysis: ${SEASON_MONTHS[@]}"

      # Extraction

      # Loop on input daily files
      IDX_NC=0
      for MAPEXT in $EXT_INS; do
          echo "Infile: $MAPEXT"
         
       # TEMPORARY LOOP TO BE REMOVED:
       if [[ -f area_tmp_${IDX_NC}.nc ]]; then
          echo "Found! I am moving to the next day.."
       else 

          # Select or compute var
          if [[ $GRID_TO_EXTRACT == "uv2t" ]]; then
             echo $(ls $MAPEXT)
             echo "Computing i=sqrt(uo*uo+vo*vo)..."
             cdo expr,'i=sqrt(uo*uo+vo*vo)' $MAPEXT name_tmp.nc
          else
             echo $(ls $MAPEXT)
             echo "cdo selname,$VAR $MAPEXT name_tmp.nc"
             cdo selname,$VAR $MAPEXT name_tmp.nc
          fi

          # Select area
          echo "cdo sellonlatbox,${LON_RANGE[0]},${LON_RANGE[1]},${LAT_RANGE[0]},${LAT_RANGE[1]} name_tmp.nc area_tmp_${IDX_NC}.nc "
          cdo sellonlatbox,${LON_RANGE[0]},${LON_RANGE[1]},${LAT_RANGE[0]},${LAT_RANGE[1]} name_tmp.nc area_tmp_${IDX_NC}.nc  
          #cdo remapnn,lat=${LAT_RANGE[0]} name_tmp.nc area_tmp_${IDX_NC}.nc
          rm name_tmp.nc       
    
       fi 

      IDX_NC=$(( $IDX_NC + 1 ))
      done
     
      # Cat extracted files
      echo "Cat extracted files.."
      cdo cat area_tmp_*.nc all_tmp.nc   
      rm area_tmp_*.nc

      # Compute the mean on the whole period
      echo "Compute the global mean.."
      cdo timmean all_tmp.nc $SECT_OUTFILE 
      
      # Compute seasonal means on the whole period
      for M_IDX in ${SEASON_MONTHS[@]}; do
           echo "Working on season: $M_IDX .."
           SECT_SEASON_OUTFILE_TMP=$( echo "$SECT_SEASON_OUTFILE_TPL" | sed -e "s/%SECT_LABEL%/$SECT_LABEL/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" -e "s/%MONTHS%/$M_IDX/g" )
           cdo selseas,${M_IDX} all_tmp.nc ${M_IDX}_tmp.nc
           cdo timmean ${M_IDX}_tmp.nc $SECT_SEASON_OUTFILE_TMP
           rm ${M_IDX}_tmp.nc
       done
       mv all_tmp.nc all_tmp.nc_${VAR}

    done

  IDX_IN=$(( $IDX_IN + 1 ))
  done

 fi


###################### POSTPROC ###########################

# Output check

# Archive


