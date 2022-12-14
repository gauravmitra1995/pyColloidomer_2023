#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

#Run all subsequent restart runs after first run using this script

simulationtype='lattice'

ntasks="1 2 3 4 5 6 7 8 9 10"
gammaA=1.0
#gammaA=0.01
gammapatch=0.0001
koninit=100.0
metropolis=1
kspring=10.0

for Np in 100;do
    for R in 50.0;do
           for rB in 1.0;do
              for Nclusters in 81;do
                 for areafraction in 0.4 0.3 0.2 0.1;do
                        for koffinit in 0.0000001 0.0001 0.01;do
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
			           			     
				     current_dir=$(pwd)
                                     
				     dir=$current_dir/$simulationtype/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${rB}/areafraction${areafraction}/epsilon${epsilon}/kspring${kspring}
				     echo $dir
				    
				     if [ -e $dir/seedlist.txt ];then
				          seedlist=$($wrapper python -u extract_seeds.py --fileprefix $dir/seedlist.txt)			     
				     fi


				     for task in $ntasks;do

                                          seed=`echo $seedlist | cut -d " " -f $task`

					  for zerodynbondlength in 'False';do

						  rundir=$($wrapper python -u update_yaml_lattice.py --Np $Np --R $R --radiusB $rB --simulationtype $simulationtype --Nclusters $Nclusters --areafraction $areafraction --seed $seed --koninit $koninit --koffinit $koffinit --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --metropolis $metropolis --zerodynbondlength $zerodynbondlength | tail -n 1)

						  echo "rundir=$rundir"

						  jobname=${simulationtype}runs_gammaA${gammaA}_gammapatch${gammapatch}_Nclus${Nclusters}_Np${Np}_R${R}_rB${rB}_areafrac${areafraction}_eps${epsilon}_seed${seed}
					   
						  cp run-all.sbatch $rundir/run-all.sbatch
 
						  cd $rundir

						  #$wrapper python -u $current_dir/run_simulation.py 					      
						  sbatch --job-name $jobname run-all.sbatch

						  cd -
			                  done
				     done
			          
			    done
                        done
                
               done 
        done
    done
done
