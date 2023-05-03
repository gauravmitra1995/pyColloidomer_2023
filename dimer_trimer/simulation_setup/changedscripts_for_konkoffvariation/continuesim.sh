#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../../dybond/run-hoomd2.9.6.bash

#Run all subsequent restart runs after first run using this script

simulationtype='polymer'
selfavoidchain=1
gammaA=0.1
metropolis=1
kspring=10.0

ntasks="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"

for Np in 100;do
    for R in 50.0;do
           for rB in 1.0;do
                for kAB in 200.0;do
                    for Nclusters in 2;do
                            for gammapatch in 0.0001 0.001 0.01;do
                                #for kon_koff in 100.0_0.0000001 50.0_0.00000005 25.0_0.000000025 20.0_0.00000002 10.0_0.00000001 5.0_0.000000005;do
                                for kon_koff in 100.0_0.001 50.0_0.0005 25.0_0.00025 20.0_0.0002 10.0_0.0001 5.0_0.00005;do
                                 for dimension in 2;do

                                    koninit=$(echo $kon_koff |cut -f 1 -d '_')
                                    koffinit=$(echo $kon_koff |cut -f 2 -d '_')

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

				     echo "Epsilon=$epsilon"
				     
                                     if [ $Nclusters -eq 1 ];then
                                          simulationtype='monomer'
                                     fi

                                     if [ $Nclusters -eq 2 ];then
                                          simulationtype='dimer'
                                     fi

                                     if [ $Nclusters -eq 3 ];then
                                          simulationtype='trimer'
                                     fi

			     
				     current_dir=$(pwd)
                                     dir=${current_dir}/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/selfavoiding/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${rB}/kAB${kAB}/epsilon${epsilon}/kon${koninit}_koff${koffinit}/dimension${dimension}/kspring${kspring}
				    
				     if [ -e $dir/seedlist.txt ];then
				          seedlist=$($wrapper python -u extract_seeds.py --fileprefix $dir/seedlist.txt)			     
				     fi
                                        

				     for task in $ntasks;do

                                          seed=`echo $seedlist | cut -d " " -f $task`
 
					  for zerodynbondlength in 'False';do

						  rundir=$($wrapper python -u update_yaml_polymer.py --Np $Np --R $R --simulationtype $simulationtype --radiusB $rB --Nclusters $Nclusters --seed $seed --koninit $koninit --koffinit $koffinit --dimension $dimension --selfavoidchain $selfavoidchain --kAB $kAB --gammaA $gammaA --gammapatch $gammapatch --metropolis $metropolis --kspring $kspring --zerodynbondlength $zerodynbondlength --dir $dir | tail -n 1)

						  echo "rundir=$rundir"
						  
						  jobname=${simulationtype}_gammaA${gammaA}_gammapatch${gammapatch}_Nclus${Nclusters}_Np${Np}_R${R}_rB${rB}_kAB${kAB}_eps${epsilon}_kon${koninit}_koff${koffinit}_kdyn${kspring}_dim${dimension}_seed${seed}
						
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
done
