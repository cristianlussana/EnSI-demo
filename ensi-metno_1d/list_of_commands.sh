#!/bin/bash

exp=$1 # 01 02 ...
i1=$2
i2=$3

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
fi

# No thinning AND No cross-validation
echo "./drive_em_all.sh -e ${exp} $str"
./drive_em_all.sh -e ${exp} $str

# Thinning AND Cross-validation
#i=$i1
#while [ "$i" -le "$i2" ]; do
#  ri=$(( ${i#0}+100 ))
#  rj=$(( ${i#0}+200 ))
#  ./drive_em_all.sh -e ${exp} ${str} -c -r ${ri} -b 20
#  for u in 25 50 75 95; do
#    ./drive_em_all.sh -e ${exp} ${str} -c -r ${ri} -b 20 -t -q ${rj} -u ${u}
#  done
#  i=$(( ${i#0}+1 ))
#done
#
exit 0
