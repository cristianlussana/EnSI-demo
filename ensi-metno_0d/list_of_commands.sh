#!/bin/bash

exp=$1 # 01 02 ...
skip=$2 # 0 1

if [ "$exp" == "01" ]; then
  str="-a 1 -s 1"
elif [ "$exp" == "02" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "03" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "04" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "05" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "06" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "07" ]; then
  str="-a 6 -s 6 -x"
elif [ "$exp" == "08" ]; then
  str="-a 3 -s 5 -x"
elif [ "$exp" == "09" ]; then
  str="-a 1 -s 6 -x"
elif [ "$exp" == "10" ]; then
  str="-a 1 -s 5 -x"
elif [ "$exp" == "11" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "12" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "13" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "14" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "15" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "16" ]; then
  str="-a 3 -s 6 -x"
elif [ "$exp" == "17" ]; then
  str="-a 3 -s 6 -x"
fi

skip_str=""
if [ "$skip" -eq "1" ]; then
  skip_str="-k"
fi

echo "./drive_em_all.sh -e $exp $str $skip_str"
./drive_em_all.sh -e $exp $str $skip_str

#
exit 0
