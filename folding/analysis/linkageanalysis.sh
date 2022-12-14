#!/bin/bash
#parameters to vary:

Np=$1
shift

R=$1
shift

radiusB=$1
shift

Nclusters=$1
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

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

metropolis=1
gammaA=0.1
gammapatch=0.0001
dt=0.0005
kspring=10.0
current_dir=$(pwd)
simulationtype='folding'
dir=$scriptdir/../simulation_setup/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${radiusB}/kAB${kAB}/epsilon${epsilon}/dimension${dimension}/kspring${kspring}


if [ -e $dir/seedlist.txt ];then
        seedlist=$($wrapper python -u extract_seeds.py --fileprefix $dir/seedlist.txt)
fi


for task in {1..300};do
   seed=`echo $seedlist | cut -d " " -f $task`

   trajectory_file=$scriptdir/../simulation_setup/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${radiusB}/kAB${kAB}/epsilon${epsilon}/dimension${dimension}/kspring${kspring}/seed${seed}/restlength${restlength}/${simulationtype}_restl${restlength}_Nc${Nclusters}_Np${Np}_R${R}_rB${radiusB}_kAB${kAB}_eps${epsilon}_dim${dimension}_kspring${kspring}_gammaA${gammaA}_gammapatch${gammapatch}_seed${seed}.allruns.gsd

   echo "trajectory_file=$trajectory_file"

   $wrapper python -u linkage_analysis.py --trajectory_file $trajectory_file --dt $dt

done


 
