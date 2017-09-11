## plot seismograms
## created by weiyw		contact: weiyw17@gmail.com
## 2017.9.10

#!/bin/bash
cd ./su_obs
# namespaces
FILE_SU=$1
############################### modify ##############################
#FILE_SU=obs_toy_vy_it1.su.shot1
#####################################################################
FILE_RSF=${FILE_SU}.rsf
sfsegyread tape=${FILE_SU} su=y endian=n format=5 > ${FILE_RSF}
LINE_REC_NUM=`sfget < ${FILE_RSF} parform=n n2 `
#awk 'BEGIN{J2=sqrt(${LINE_REC_NUM});}'
J2=` echo "sqrt(${LINE_REC_NUM}) - 1" | bc `
echo ${J2}
sfwindow j2=${J2} < ${FILE_RSF} | sfgrey plicp=97 | sfpen
