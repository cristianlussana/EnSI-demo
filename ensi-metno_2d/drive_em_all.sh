#!/bin/bash
#-----------------------------------------------------------------

# Include here the instructions needed to properly run R on your machine
source /home/klinogrid/etc/profile_rhel8_r413.sh

# Get input parameters
cv=0
thin=0
n_adj=0
cv_perc=9999
cv_setseed=9999
thin_perc=9999
thin_setseed=9999
exp=exp00
dir_main=./
ensi=$dir_main/ensi-metno_2d.r
dir_in=./data
extra=0
yextra=0
while getopts "s:a:e:xycb:r:tq:u:" Option
do
  case $Option in
  s ) n_steps=$OPTARG
  ;;
  a ) n_adj=$OPTARG
  ;;
  e ) exp=$OPTARG
  ;;
  x ) extra=1
  ;;
  y ) yextra=1
  ;;
  c ) cv=1
  ;;
  r ) cv_setseed=$OPTARG
  ;;
  b ) cv_perc=$OPTARG
  ;;
  t ) thin=1
  ;;
  q ) thin_setseed=$OPTARG
  ;;
  u ) thin_perc=$OPTARG
  ;;
  esac
done

# Elaboration
if [ "$cv" -eq "0" ] && [ "$thin" -eq "0" ]; then
  res_str="res_nocv_nothin"
  tmp_str="tmp_nocv_nothin"
  cv_str=""
  thin_str=""
elif [ "$cv" -eq "0" ] && [ "$thin" -ne "0" ]; then
  res_str="res_nocv_thin${thin_perc}_${thin_setseed}"
  tmp_str="tmp_nocv_thin${thin_perc}_${thin_setseed}"
  cv_str=""
  thin_str="--thinobs --thinobs_setseed $thin_setseed --thinobs_perc $thin_perc"
elif [ "$cv" -ne "0" ] && [ "$thin" -eq "0" ]; then
  res_str="res_cv${cv_perc}_${cv_setseed}_nothin"
  tmp_str="tmp_cv${cv_perc}_${cv_setseed}_nothin"
  cv_str="--cvobs --cvobs_setseed $cv_setseed --cvobs_perc $cv_perc"
  thin_str=""
elif [ "$cv" -ne "0" ] && [ "$thin" -ne "0" ]; then
  res_str="res_cv${cv_perc}_${cv_setseed}_thin${thin_perc}_${thin_setseed}"
  tmp_str="tmp_cv${cv_perc}_${cv_setseed}_thin${thin_perc}_${thin_setseed}"
  cv_str="--cvobs --cvobs_setseed $cv_setseed --cvobs_perc $cv_perc"
  thin_str="--thinobs --thinobs_setseed $thin_setseed --thinobs_perc $thin_perc"
fi
dir_tmp=./output/ensi-metno_2d/case_20230807/exp${exp}/$tmp_str
dir_res=./output/ensi-metno_2d/case_20230807/exp${exp}/$res_str
dir_past=$dir_tmp/past
dir_fresh=$dir_tmp/fresh
ffin_l_template=$dir_fresh/analysis_%V_2d_%Y%m%dT%HZ.rdata
ffin_r_template=$dir_past/analysis_%V_2d_%Y%m%dT%HZ.rdata
ffout_template=$dir_fresh/analysis_%V_2d_%Y%m%dT%HZ.rdata
[[ -d $dir_res ]] && rm -r $dir_res
[[ -d $dir_past ]] && rm -r $dir_past
[[ -d $dir_fresh ]] && rm -r $dir_fresh
mkdir -p $dir_res
mkdir -p $dir_past
mkdir -p $dir_fresh
for hh in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do
  cp $dir_in/MEPSpp_RR1-MEPS_lcc_2d_20230807T${hh}Z.rdata $dir_past/analysis_RR1_2d_20230807T${hh}Z.rdata
  cp $dir_in/MEPSpp_TA-MEPS_lcc_2d_20230807T${hh}Z.rdata $dir_past/analysis_TA_2d_20230807T${hh}Z.rdata
  cp $dir_in/MEPSpp_RR1-MEPS_lcc_2d_20230807T${hh}Z.rdata $dir_fresh/analysis_RR1_2d_20230807T${hh}Z.rdata
  cp $dir_in/MEPSpp_TA-MEPS_lcc_2d_20230807T${hh}Z.rdata $dir_fresh/analysis_TA_2d_20230807T${hh}Z.rdata
done
datel1="2023-08-07"

if [ "$yextra" -eq "1" ]; then
  adjcont=1
  stepcont=1
  for hh in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do
    datel=${datel1}T${hh}
    for varl in RR1 TA; do
      echo "+-----------------------------------------------------------+"
      echo "$ensi --config_file ${dir_main}/exp${exp}/exp${exp}_extra.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template"
      $ensi --config_file ${dir_main}/exp${exp}/exp${exp}_extra.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_Str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template
    done
  done
  cp $dir_fresh/*.rdata $dir_past/
fi

adjcont=1
while [ "$adjcont" -le "$n_adj" ]; do 
  for hh in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do
    datel=${datel1}T${hh}
    for varl in RR1 TA; do
      stepcont=1
      while [ "$stepcont" -le "$n_steps" ]; do
        echo "+-----------------------------------------------------------+"
        echo "$ensi --config_file ${dir_main}/exp${exp}/exp${exp}.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template"
        $ensi --config_file ${dir_main}/exp${exp}/exp${exp}.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template
        stepcont=$(( ${stepcont#0}+1 ))
      done
    done
  done
  cp $dir_fresh/*.rdata $dir_past/
  adjcont=$(( ${adjcont#0}+1 ))
done

if [ "$extra" -eq "1" ]; then
  adjcont=1
  stepcont=1
  for hh in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do
    datel=${datel1}T${hh}
    for varl in RR1 TA; do
      echo "+-----------------------------------------------------------+"
      echo "$ensi --config_file ${dir_main}/exp${exp}/exp${exp}_extra.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template"
      $ensi --config_file ${dir_main}/exp${exp}/exp${exp}_extra.ini --varl $varl --datel $datel --adjust $adjcont --step $stepcont $cv_str $thin_Str --ffin_l_template $ffin_l_template --ffin_r_template $ffin_r_template --ffout_template $ffout_template
    done
  done
fi

# Save results
cp -r $dir_fresh/*.rdata $dir_res/
rm -r $dir_tmp

#--------------------------
# Clean and Exit 
#--------------------------
duration=$SECONDS
echo "the end: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
exit 0
