#!/bin/bash
#parameters to vary:

restlength=$1
shift

R=$1
shift

radiusB=$1
shift

Np=$1
shift

areafraction=$1
shift

koninit=$1
shift

koffinit=$1
shift

gammaA=$1
shift

gammapatch=$1
shift


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

Nclusters=81
kspring=10.0
metropolis=1

module purge
scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

#This combines data from all seeds to obtain the histograms and the plots of valence/structure over time
$wrapper python -u averageoverallseeds_valenceandstructures.py --Nclusters $Nclusters --R $R --radiusB $radiusB --Np $Np --areafraction $areafraction --epsilon $epsilon --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --restlength $restlength 

#This averages data over all seeds and finds the SD
$wrapper python -u averageoverallseeds_valenceandstructures_geterrorbarsforseeds.py --Nclusters $Nclusters --R $R --radiusB $radiusB --Np $Np --areafraction $areafraction --epsilon $epsilon --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --restlength $restlength

$wrapper python -u histogram_valences_overallseeds.py --Nclusters $Nclusters --R $R --radiusB $radiusB --Np $Np --areafraction $areafraction --epsilon $epsilon --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --restlength $restlength

$wrapper python -u histogram_structures_overallseeds.py --Nclusters $Nclusters --R $R --radiusB $radiusB --Np $Np --areafraction $areafraction --epsilon $epsilon --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --restlength $restlength

$wrapper python -u averageoverallseeds_linearchainlengths.py --Nclusters $Nclusters --R $R --radiusB $radiusB --Np $Np --areafraction $areafraction --epsilon $epsilon --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --restlength $restlength

