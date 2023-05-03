#!/bin/bash
#parameters to vary:

Np=$1
shift

R=$1
shift

radiusB=$1
shift

gammapatch=$1
shift

kAB=$1
shift

dimension=$1
shift

koninit=$1
shift

koffinit=$1
shift

restlength=$1
shift

Nclusters=2

if [ $koffinit -eq 0 ];then
       epsilon='infinite'
else
      if [ $koninit -eq 0 ];then
             epsilon=0.0
      else
             epsilon_raw=$(echo "l($koninit/$koffinit)" | bc -l)
             epsilon=`printf "%.1f" $epsilon_raw`
      fi
fi

if [ $Nclusters -eq 2 ];then
        simulationtype='dimer'
fi

if [ $Nclusters -eq 3 ];then
        simulationtype='trimer'
fi


metropolis=1
gammaA=0.1
dt=0.001
kspring=10.0
current_dir=$(pwd)

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../../dybond/run-hoomd2.9.6.bash

$wrapper python -u averageoverallseeds.py --Nclusters $Nclusters --Np $Np --R $R --radiusB $radiusB --kAB $kAB --dimension $dimension --epsilon $epsilon --kspring $kspring --gammaA $gammaA --gammapatch $gammapatch --r0 $restlength --koninit $koninit --koffinit $koffinit


 
