#!/bin/bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1

#SBATCH --mem 4GB

#SBATCH -t 12:00:00
##SBATCH --dependency=singleton

##SBATCH --gres=gpu:1
##SBATCH --gres=gpu:1g.10gb:1


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

module purge
scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

metropolis=1
dt=0.001
kspring=10.0
Nclusters=81

current_dir=$(pwd)
simulationtype='lattice'
dir=$scriptdir/../simulation_setup/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${radiusB}/areafraction${areafraction}/epsilon${epsilon}/kspring${kspring}
echo $dir

if [ -e $dir/seedlist.txt ];then
        seedlist=`cat $dir/seedlist.txt`
fi

echo $seedlist 

ntasks="1 2 3 4 5 6 7 8 9 10"

for task in $ntasks;do

   seed=`echo $seedlist | cut -d " " -f $task`

   fileprefix=$scriptdir/../simulation_setup/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${radiusB}/areafraction${areafraction}/epsilon${epsilon}/kspring${kspring}/seed${seed}/restlength${restlength}/${simulationtype}_restl${restlength}_Nc${Nclusters}_Np${Np}_R${R}_rB${radiusB}_areafrac${areafraction}_eps${epsilon}_kspring${kspring}_gammaA${gammaA}_gammapatch${gammapatch}_seed${seed}
   $wrapper python -u gsdfile_forsingleframe.py --fileprefix $fileprefix


   trajectory_file=$scriptdir/../simulation_setup/finalframes/${simulationtype}_restl${restlength}_Nc${Nclusters}_Np${Np}_R${R}_rB${radiusB}_areafrac${areafraction}_eps${epsilon}_kspring${kspring}_gammaA${gammaA}_gammapatch${gammapatch}_seed${seed}.finalframe.gsd

   echo "trajectory_file=$trajectory_file"

   $wrapper python -u unwrapping_trajectories.py --trajectory_file $trajectory_file --dt $dt

done



