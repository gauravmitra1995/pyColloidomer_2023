#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

#Run all subsequent restart runs after first run using this script 

simulationtype='folding'
gammaA=0.1
gammapatch=0.0001
koninit=200.0
metropolis=1
kspring=10.0

kT=1.3 #run alternatively at high temp and then again at low temp kT=1.0, continue this cycle for as many times as you want. Just set the right temperature every time before submitting a restart run. 

for Np in 200;do
    for R in 20.0;do
           for rB in 1.0;do
                for kAB in 200.0;do
                    for Nclusters in 7;do
                            for koffinit in 2.0;do
                               for dimension in 2;do
				       
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
				     dir=${current_dir}/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${rB}/kAB${kAB}/epsilon${epsilon}/dimension${dimension}/kspring${kspring}
				     #seedlist="12197099"
				     if [ -e $dir/seedlist.txt ];then
				          seedlist=$($wrapper python -u extract_seeds.py --fileprefix $dir/seedlist.txt)			     
				     fi

				     echo $seedlist
                                     
                                     #for task in {1..1};do	
				     for task in {1..300};do
                                               seed=`echo $seedlist | cut -d " " -f $task`

                                               for zerodynbondlength in 'False';do

                                                 rundir=$($wrapper python -u update_yaml_folding.py --Np $Np --R $R --simulationtype $simulationtype --radiusB $rB --Nclusters $Nclusters --seed $seed --koninit $koninit --koffinit $koffinit --dimension $dimension --kAB $kAB --gammaA $gammaA --gammapatch $gammapatch --metropolis $metropolis --kspring $kspring --zerodynbondlength $zerodynbondlength --kT $kT | tail -n 1)

                                                 echo "rundir=$rundir"

						 jobname=folding_gammaA${gammaA}_gammapatch${gammapatch}_Nclus${Nclusters}_Np${Np}_R${R}_rB${rB}_kAB${kAB}_eps${epsilon}_kdyn${kspring}_dim${dimension}_seed${seed}

                                                 cp run-all.sbatch $rundir/run-all.sbatch

                                                 cd $rundir

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
done
