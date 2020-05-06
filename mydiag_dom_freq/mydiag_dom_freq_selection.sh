# !/bin/bash
# 
#set -x
set -e 
#set -u
#
# by
# AC Goglio (annachiara.goglio@cmcc.it)
# CMCC Bologna 
# 08/04/2020
#
# HOW TO RUN THIS SCRIPT: Set the variables in the mydiag_dom_freq.ini file reading the instructions and just run it!
# 
# ==========================================
# Reading ini file
# ==========================================
echo "Reading the ini file.."
ls 
INIFILE="mydiag_dom_freq.ini"
if [[ -f ./${INIFILE} ]]; then
   source ${INIFILE}
   echo "..Done!"
else
   echo "I cannot find the ini file ${INIFILE}..Why?! The file MUST be stored in the same dir of this script.. "
fi
#

###########################################
###########################################
# DO NOT CHANGE THE CODE BENEATH THIS LINE
# UNLESS YOU KNOW WHAT YOU ARE DOING!!
###########################################
###########################################

# Exp infos check
if [[ -d ${EXP_PATH} ]] && [[ -d ${WORK_DIR}/${EXP_NAME} ]] ; then
   echo "I am working on exp: ${EXP_NAME}"
   echo "the src directory of the exp is: ${EXP_PATH}"
   echo "the work dir of the exp is: ${WORK_DIR}"
   sleep 5
elif [[ ! -d ${EXP_PATH} ]]; then
   echo "ERROR: NOT found ${EXP_PATH}... WHY? Check EXP_PATH and EXP_NAME and try again! "
   exit
elif [[ ! -d ${WORK_DIR}/${EXP_NAME} ]]; then
   echo "ERROR: exp ${EXP_NAME} NOT found in ${WORK_DIR}... WHY? Check WORK_DIR and EXP_NAME and try again! "
   exit
fi

# Domain check and setting
if [[ $DOMAIN == 'med' ]]; then
   DOMAIN_NAME=$DOMAIN
   DOMAIN=" [  [[-5.99 , 1 ], [34,42]] , [[0,36.3],[30.19,45.99]] ] "
   echo "Working on Mediterranean Basin.."
   echo "Coordinates: $DOMAIN"
   sleep 3
elif [[ $DOMAIN == 'atl' ]]; then
   DOMAIN_NAME=$DOMAIN
   DOMAIN="[ [[-5.99 , 1 ], [42.5,45.99]] , [[-18.12,-5.99],[30.19,45.99]] ]"
   echo "Working on Atlantic Box.."
   echo "Coordinates: $DOMAIN"
   sleep 3
else
   DOMAIN_NAME="DomainUserDefn"
   echo "Working on a user defined domain: good luck!"
   echo "Coordinates: $DOMAIN"
   sleep 3   
fi

# Frequency warnings and parameters setting
if [[ $FREQ == '1d' ]]; then
   echo "Working on DAILY model outputs"
   #
   T_TEMPLATEFILE="*_1d_????????_grid_T.nc"
   U_TEMPLATEFILE="*_1d_????????_grid_U.nc"
   V_TEMPLATEFILE="*_1d_????????_grid_V.nc"
   T_TEMPLATEFILE4MASK=${T_TEMPLATEFILE}
   #
   KEY_VOTEMPER="tem_t_all_l"
   KEY_VOSALINE="tem_s_all_l"
   KEY_2D="tem_t"
   #
   REDUCE_COMMAND_LINE='bsub -J `basename $0` -o log_out.%J -e log_err.%J -q serial_6h -R "rusage[mem=5000]" "$FilterFile $ExpDir/exp_mrsp.sh/`basename $XmlFile` | mapreduce.py $XmlFile"'
   REDUCED_LN_LIST='make_list $wdir tra_t_gb_ts.nc ; make_list $wdir tra_p_gb_ts.nc ; make_list $wdir tra_n_gb_ts.nc ; make_list $wdir tra_t_sc_ts.nc ; make_list $wdir tra_p_sc_ts.nc ; make_list $wdir tra_n_sc_ts.nc ; make_list $wdir tra_t_ot_ts.nc ; make_list $wdir tra_p_ot_ts.nc ; make_list $wdir tra_n_ot_ts.nc ; make_list $wdir tra_t_co_ts.nc ; make_list $wdir tra_p_co_ts.nc ; make_list $wdir tra_n_co_ts.nc ; make_list $wdir tra_t_me_ts.nc ; make_list $wdir tra_p_me_ts.nc ; make_list $wdir tra_n_me_ts.nc ; make_list $wdir T.nc_SST_ts.nc ; make_list $wdir T.nc_0_150_ts.nc ; make_list $wdir T.nc_150_600_ts.nc ; make_list $wdir T.nc_600_btm_ts.nc ; make_list $wdir T.nc_basin_ts.nc ; make_list $wdir S.nc_SSS_ts.nc ; make_list $wdir S.nc_0_150_ts.nc ; make_list $wdir S.nc_150_600_ts.nc ; make_list $wdir S.nc_600_btm_ts.nc ; make_list $wdir S.nc_basin_ts.nc ; make_list $wdir sossheig_ts.nc ; make_list $wdir somxl010_ts.nc ; make_list $wdir sohefldo_ts.nc ; make_list $wdir sowaflup_ts.nc ; make_list $wdir soevapor_ts.nc ; make_list $wdir soprecip_ts.nc ; make_list $wdir sorunoff_ts.nc ; make_list $wdir soshfldo_ts.nc ; make_list $wdir solofldo_ts.nc ; make_list $wdir sosefldo_ts.nc ; make_list $wdir solafldo_ts.nc ; make_list $wdir velmodw_ts.nc ; make_list $wdir V.nc_SSV_ts.nc ; make_list $wdir V.nc_0_150_ts.nc ; make_list $wdir V.nc_150_600_ts.nc ; make_list $wdir V.nc_600_btm_ts.nc ; make_list $wdir V.nc_basin_ts.nc ; make_list $wdir K.nc_SSK_ts.nc ; make_list $wdir K.nc_0_150_ts.nc ; make_list $wdir K.nc_150_600_ts.nc ; make_list $wdir K.nc_600_btm_ts.nc ; make_list $wdir K.nc_basin_ts.nc'
   #
   if [[ ${DOMAIN_NAME} == "med" ]] || [[ ${DOMAIN_NAME} == "DomainUserDefn" ]] ; then
      VOZOCRTX="%t_cur_u_l;%t_sec_u_all,1" 
      KEY_VOZOCRTX="cur_u"
      VOMECRTY="%t_cur_v_l;%t_sec_v_all,1"
      KEY_VOMECRTY="cur_v"
      SOZOTAUX="%t_cur_uw" 
      KEY_SOZOTAUX="cur_u"
      SOMETAUY="%t_cur_vw"
      KEY_SOMETAUY="cur_v"
   else
      VOZOCRTX="%t_cur_u_l"
      KEY_VOZOCRTX="cur_u_ko"
      VOMECRTY="%t_cur_v_l"
      KEY_VOMECRTY="cur_v_ko"
      SOZOTAUX="%t_cur_uw"
      KEY_SOZOTAUX="cur_u_ko"
      SOMETAUY="%t_cur_vw"
      KEY_SOMETAUY="cur_v_ko"
   fi
elif [[ $FREQ == '1h' ]]; then
   echo "Working on HOURLY model outputs"
   echo "WARNING: at least a daily output file is needed.. "
   #
   T_TEMPLATEFILE="*_1h_????????_grid_T.nc"
   U_TEMPLATEFILE="*_1h_????????_grid_U.nc"
   V_TEMPLATEFILE="*_1h_????????_grid_V.nc"
   T_TEMPLATEFILE4MASK="*_1d_????????_grid_T.nc"
   #
   REDUCE_COMMAND_LINE='bsub -J `basename $0` -o log_out.%J -e log_err.%J -q serial_6h -R "rusage[mem=5000]" "$FilterFile $ExpDir/exp_mrsp.sh/`basename $XmlFile` | grep "grid_T" | mapreduce.py $XmlFile"'
   REDUCED_LN_LIST='make_list $wdir sossheig_ts.nc'
   # 1h fields are 2d so all the 3d fileds analysis must be avoided:
   KEY_VOTEMPER="tem_t_all_l_ko"
   KEY_VOSALINE="tem_s_all_l_ko"
   KEY_2D="tem_t_ko"
   #
   VOZOCRTX="%t_cur_u_l"
   KEY_VOZOCRTX="cur_u_ko"
   VOMECRTY="%t_cur_v_l"
   KEY_VOMECRTY="cur_v_ko"
   SOZOTAUX="%t_cur_uw"
   KEY_SOZOTAUX="cur_u_ko"
   SOMETAUY="%t_cur_vw"
   KEY_SOMETAUY="cur_v_ko"
else
   echo "Check FREQ value and try again!"
   exit
fi

# Vertical transect settings
if [[ $TRA_V_SECTION == "GB" ]]; then
   # If you want to change the depth of the section change the last two number..
   V_TRANSECT_TRA="[  [[-5.48,-5.47 ], [35,37]] , ["d", [[0.5,700]] ] ]"
   echo "I am computing the Gibraltar transport (Med side).."
elif [[ $TRA_V_SECTION == "GBA" ]]; then
   # If you want to change the depth of the section change the last two number..
   V_TRANSECT_TRA="[  [[-5.48,-5.47 ], [35,37]] , ["d", [[0.5,700]] ] ]"
   echo "I am computing the Gibraltar transport (Atlantic side).."
elif [[ $TRA_V_SECTION == "DARD" ]]; then
   # If you want to change the depth of the section change the last two number..
   V_TRANSECT_TRA="[  [[26.14,26.16 ], [39.97,40.11]] , ["d", [[0.5,600]] ] ]"
   echo "I am computing the Dardanelles transport.."
fi

# Clean old runs 
TO_BE_REMOVED=("sossheig.nc" "votemper.nc" "/pack/med_prod/diag_base.xml" )
echo "I am going to REMOVE the following temporary files from exp directory that must be rebuilt: ${TO_BE_REMOVED}"
sleep 2
for FILE_2_RM in ${TO_BE_REMOVED}; do
    if [[ -f ${EXP_PATH}/${FILE_2_RM} ]]; then
       rm -v ${EXP_PATH}/${FILE_2_RM} 
    fi
done

# Built the scripts from the templates with respect to the choosen settings  
echo "I am going to BUILD the required scripts into the EXP_PATH/EXP_NAME/pack directory.. In order to avoid the removal of old versions, if found, this will be moved to *_old files.."
TO_BE_BUILT=("diag_base_eas3.TEMPLATE.xml" "mydiag.sh" "gen_xml.sh" "exp_mrsp.sh" "exp_filter_reduce.sh" "exp_filter_map.sh" )
echo "${TO_BE_BUILT[@]}"

# Loop on single files to be built
for FILE_2_BUILT in ${TO_BE_BUILT[@]}; do
    echo "I am working on file ${FILE_2_BUILT}.."
    
    # Where to put new files
    if [[ ${FILE_2_BUILT} == "diag_base_eas3.TEMPLATE.xml" ]]; then
       WHERE_2_BUILT="pack/med_prod"
    else
       WHERE_2_BUILT="pack/med_prod/bin"
    fi
    echo "The file is going to be created in ${EXP_PATH}/${WHERE_2_BUILT}/"
    
    # Mv the old versions of files (just in case..)
    if [[ -f ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT} ]]; then
       echo "An old version of this file has been found.. I am moving it to _old file.. (Just in case!!)"
       #mv -v ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT} ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT}_old
       #rm -v ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT} ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT}
       echo "Done!"
    fi

    # Sed file creation and sobstitution of parematers in the templates  
    SED_FILE=sed_file.txt
    cat << EOF > ./${SED_FILE}
   s/%DOMAIN%/${DOMAIN//\//\\/}/g
   #
   s/%T_TEMPLATEFILE%/${T_TEMPLATEFILE//\//\\/}/g
   s/%U_TEMPLATEFILE%/${U_TEMPLATEFILE//\//\\/}/g
   s/%V_TEMPLATEFILE%/${V_TEMPLATEFILE//\//\\/}/g
   s/%T_TEMPLATEFILE4MASK%/${T_TEMPLATEFILE4MASK//\//\\/}/g
   #
   s/%V_TRANSECT_TRA%/${V_TRANSECT_TRA//\//\\/}/g
   #
   s/%VOZOCRTX%/${VOZOCRTX//\//\\/}/g
   s/%KEY_VOZOCRTX%/${KEY_VOZOCRTX//\//\\/}/g
   s/%VOMECRTY%/${VOMECRTY//\//\\/}/g
   s/%KEY_VOMECRTY%/${KEY_VOMECRTY//\//\\/}/g
   s/%SOZOTAUX%/${SOZOTAUX//\//\\/}/g
   s/%KEY_SOZOTAUX%/${KEY_SOZOTAUX//\//\\/}/g   
   s/%SOMETAUY%/${SOMETAUY//\//\\/}/g   
   s/%KEY_SOMETAUY%/${KEY_SOMETAUY//\//\\/}/g
   s/%KEY_VOTEMPER%/${KEY_VOTEMPER//\//\\/}/g
   s/%KEY_VOSALINE%/${KEY_VOSALINE//\//\\/}/g
   s/%KEY_2D%/${KEY_2D//\//\\/}/g
   #
   s/%REDUCE_COMMAND_LINE%/${REDUCE_COMMAND_LINE//\//\\/}/g
   s/%REDUCED_LN_LIST%/${REDUCED_LN_LIST//\//\\/}/g
EOF

      sed -f ${SED_FILE} ./templates/${FILE_2_BUILT}_template > ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT}
      rm ${SED_FILE}

      # Chenge user permissions
      chmod u+x ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT}
done

# Final check on mydiag settings:
# Loop on single scripts that must have been built
for FILE_2_BUILT in ${TO_BE_BUILT[@]}; do
    echo "I am checking on file ${FILE_2_BUILT}.."

    # Where to look for the new files
    if [[ ${FILE_2_BUILT} == "diag_base_eas3.TEMPLATE.xml" ]]; then
       WHERE_2_BUILT="pack/med_prod"
    else
       WHERE_2_BUILT="pack/med_prod/bin"
    fi

    # file existence check
    if [[ -f ${EXP_PATH}/${WHERE_2_BUILT}/${FILE_2_BUILT} ]]; then
       echo "OK: this file has been found!"
    else
       echo "ERROR: file NOT found!!! Why?"
       echo "mydiag cannot be executed without it: I am exiting.."
       exit
    fi
done
##################################################
# MYDIAG DIAGNOSTIC PACK RUNNING
##################################################
# mydiag.sh automatically launched by this script
if [[ ${MYDIAG_FLAG} == 1 ]]; then
 # Moving to the exp src directory and loading the environment..
 echo "I am moving to the exp src directory: ${EXP_PATH}"
 cd ${EXP_PATH}
 echo "I am loading the environment.."
 . setup.sh
 echo "..Done!"
 # Running mydiag script
 echo "Running mydiag.sh script: good luck! (Check the queue.. qstat )"
 mydiag.sh ${WORK_DIR}/${EXP_NAME}
 sleep 2
 qstat
 
  if [[ $FREQ == '1h' ]]; then 
    echo "You choose to work on hourly model outputs. For some reasons in this case the code does not cat the montly time-series in a single file."
    echo "It in thus suggested to do the following steps when the script has finished:"
    echo "module load CDO/cdo-1.5.9 (WARNING: pay attention to the cdo version you are loading, successive ones do not work on this files..)"
    echo "cdo mergetime ${WORK_DIR}/${EXP_NAME}/exp_mrsp.sh/diag_base.xml/20????/sossheig_ts.nc ${WORK_DIR}/${EXP_NAME}/exp_mrsp.sh/diag_base.xml/sossheig_ts.nc "
  elif [[ $FREQ == '1d' ]]; then
    echo "For some reasons in this case the code needs a second run in order to cat the montly time-series in a single file."
    echo "It in thus suggested to do the following steps when the script has finished:"
    echo "mydiag.sh ${WORK_DIR}/${EXP_NAME}"
  fi 
# mydiag.sh handled by the user
else
 echo "DONE: now you can run mydiag.sh on your exp for src directory"
 echo ". setup.sh"
 echo "mydiag.sh ${WORK_DIR}/${EXP_NAME}"
fi
#
echo "Thanks for trusting this script!"
echo "Bye"
