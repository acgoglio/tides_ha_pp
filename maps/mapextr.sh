#!/bin/bash
#
# ACG 04/11/2019
# Script for TS extraction  
# Ini file: map_extr.ini 
#
#set -u
set -e
#set -x 
################### PREPROC ###################################

# Source ini file
  INIFILE='map_extr.ini'
  source $INIFILE

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
     cp ${SRC_DIR}/$INIFILE ${ANA_WORKDIR}

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
        if [[ $GRID_TO_EXTRACT != "uv2t" ]]; then
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
  
  if [[ $MAP_FLAG == 1 ]]; then
    echo "Map analisys required!"
    module load $MAP_MODULE
    echo "Loading the Enviroment.."
  fi


####################### ANALYSIS ##############################

# Extract values from netCDF

 if [[ $MAP_FLAG == 1 ]]; then

  echo "------ Maps Extraction from model outputs ------"

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
     
     ## Lop on vertical levels
     #VL_IDX=0
     #echo "# of v-lev = ${#VAR3D_LEV[@]}"
     #while [[ $VL_IDX -lt ${#VAR3D_LEV[@]} ]]; do

       #VLEV=${VAR3D_LEV[$VL_IDX]}
       VLEV="allv"
       #echo "Depth: $VLEV .. "

       MAP3D_OUTFILE=$( echo "$MAP3D_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

     #  for M_IDX in ${SEASON_MONTHS[@]}; do
     #      STRING_TMP=$( echo "$MAP3D_SEASON_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" -e "s/%MONTHS%/$M_IDX/g" )
     #      echo ${STRING_TMP}
     #      declare MAP3D_OUTFILE_${M_IDX}="$( echo ${STRING_TMP} )"
     #  done

      echo "# MAP extraction from ${ANA_INTAG[$IDX_IN]}" 
      echo "# VAR: ${VAR}" 
      echo "# LEV: ${VLEV}" 
      echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"
      echo "# Seasonal analysis: ${SEASON_MONTHS[@]}"

      # Extraction

      # Loop on input daily files
      IDX_NC=0
      for MAPEXT in $EXT_INS; do
          echo "Infile: $MAPEXT"
          
          # Select or compute var
          if [[ $GRID_TO_EXTRACT == "uv2t" ]]; then
             echo "Computing i=sqrt(uo*uo+vo*vo)..."
             #cdo expr,'i=sqrt(uo*uo+vo*vo)' $MAPEXT name_tmp_${IDX_NC}.nc
             cdo expr,'ut=uo;vt=vo' $MAPEXT name_tmp_${IDX_NC}.nc
          else
             echo "Selecting VAR:$VAR..." 
             cdo selname,$VAR $MAPEXT name_tmp_${IDX_NC}.nc
          fi

          # Select vertical level
          #cdo intlevel,${VLEV} name_tmp.nc lev_tmp_${IDX_NC}.nc  
          #rm name_tmp.nc       
 

      IDX_NC=$(( $IDX_NC + 1 ))
      done
     
      # Cat extracted files
      echo "Cat extracted files.."
      cdo cat name_tmp_*.nc all_tmp.nc   
      rm name_tmp_*.nc

      # Compute the mean on the whole period
      echo "Compute the global mean.."
      cdo timmean all_tmp.nc $MAP3D_OUTFILE 

      # Compute seasonal means on the whole period
      for M_IDX in ${SEASON_MONTHS[@]}; do
           echo "Working on season: $M_IDX .."
           MAP3D_SEASON_OUTFILE_TMP=$( echo "$MAP3D_SEASON_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" -e "s/%MONTHS%/$M_IDX/g" )
           cdo selseas,${M_IDX} all_tmp.nc ${M_IDX}_tmp.nc
           cdo timmean ${M_IDX}_tmp.nc $MAP3D_SEASON_OUTFILE_TMP
           #rm all_tmp.nc
           rm ${M_IDX}_tmp.nc
       done

     # for SEASON in ${SEASON_MONTHS[@]}; do
     #     cdo selseas,$SEASON all_tmp.nc ${SEASON}_tmp.nc
     #     echo "timmean ${SEASON}_tmp.nc MAP3D_OUTFILE_${SEASON}"
     #     cdo timmean ${SEASON}_tmp.nc ${MAP3D_OUTFILE_${SEASON}}
     #     rm all_tmp.nc
     #     rm ${SEASON}_tmp.nc
     # done



     #VL_IDX=$(( $VL_IDX + 1 ))
     #done

    rm all_tmp.nc
    done

    ############# 2D VARS #########################

    echo "I am working on 2D vars..."

    # Loop on 2D FIELDS
    for VAR in ${VAR2D_NAME[@]}; do
     echo "Extracting field: $VAR .. "

     # Input file names (linked in work dir)
     EXT_INS=$(  echo "${ANA_INFILES_TPL[${IDX_IN}]}" | sed -e "s/%YYYYMMDD%/"*"/g" )

       VLEV=0
       echo "Depth: $VLEV .. "

       MAP2D_OUTFILE=$( echo "$MAP2D_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" )

      echo "# MAP extraction from ${ANA_INTAG[$IDX_IN]}" 
      echo "# VAR: ${VAR}" 
      echo "# LEV: ${VLEV}" 
      echo "# Period: ${ANA_STARTDATE}-${ANA_ENDDATE}"
      echo "# Seasonal analysis: ${SEASON_MONTHS[@]}"

      # Extraction

      # Loop on input daily files
      IDX_NC=0
      for MAPEXT in $EXT_INS; do
          echo "Infile: $MAPEXT"

          # Select var
          cdo selname,$VAR $MAPEXT name_tmp_${IDX_NC}.nc

      IDX_NC=$(( $IDX_NC + 1 ))
      done

      # Cat extracted files
      echo "Cat extracted files.."
      cdo cat name_tmp_*.nc all_tmp.nc
      rm name_tmp_*.nc

      # Compute the mean on the whole period
      echo "Compute the global mean.."
      cdo timmean all_tmp.nc $MAP2D_OUTFILE

      # Compute seasonal means on the whole period

       for M_IDX in ${SEASON_MONTHS[@]}; do
           echo "Working on season: $M_IDX .."
           MAP2D_SEASON_OUTFILE_TMP=$( echo "$MAP2D_SEASON_OUTFILE_TPL" | sed -e "s/%LEV%/$VLEV/g" -e "s/%FIELD%/${VAR}/g" -e "s/%ANA_STARTDATE%/${ANA_STARTDATE}/g" -e "s/%ANA_ENDDATE%/${ANA_ENDDATE}/g" -e "s/%INDATASET%/${ANA_INTAG[$IDX_IN]}/g" -e "s/%MONTHS%/$M_IDX/g" )
           cdo selseas,$M_IDX all_tmp.nc ${M_IDX}_tmp.nc
           cdo timmean ${M_IDX}_tmp.nc $MAP2D_SEASON_OUTFILE_TMP
           rm ${M_IDX}_tmp.nc
       done

    rm all_tmp.nc
    done

   IDX_IN=$(( $IDX_IN + 1 ))
   done

 fi


###################### POSTPROC ###########################

# Output check

# Archive


