#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

ntasks=300

#Very first run set up using this script

simulationtype='folding'
gammaA=0.1
gammapatch=0.0001
koninit=200.0
metropolis=1
kspring=10.0
kT=1.0  

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
				  
	                          seedlist=""
				  for ((x=1; x<=$ntasks; x++))
				  do
					  no=`od -vAn -N3 -tu < /dev/urandom`
					  seedlist=`echo "$seedlist ${no}"`
				  done


				  dir=${current_dir}/${simulationtype}/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${rB}/kAB${kAB}/epsilon${epsilon}/dimension${dimension}/kspring${kspring}
				  mkdir -p $dir
                                  

				  if [ -e $dir/seedlist.txt ];then
					  echo -n " " $seedlist >> $dir/seedlist.txt
				  else
					  echo -n $seedlist > $dir/seedlist.txt
				  fi

                                  for seed in $seedlist;do

					       for zerodynbondlength in 'False';do

						 rundir=$($wrapper python -u update_yaml_folding.py --Np $Np --R $R --simulationtype $simulationtype --radiusB $rB --Nclusters $Nclusters --seed $seed --koninit $koninit --koffinit $koffinit --dimension $dimension --kAB $kAB --gammaA $gammaA --gammapatch $gammapatch --metropolis $metropolis --kspring $kspring --zerodynbondlength $zerodynbondlength --kT $kT | tail -n 1)


						 echo "rundir=$rundir"

						 jobname=folding_gammaA${gammaA}_gammapatch${gammapatch}_Nclus${Nclusters}_Np${Np}_R${R}_rB${rB}_kAB${kAB}_eps${epsilon}_kdyn${kspring}_dim${dimension}_seed${seed}

						 cp run-all.sbatch $rundir/run-all.sbatch
					
						 cd $rundir

						 $wrapper python -u $current_dir/run_simulation.py
						 #sbatch --job-name $jobname run-all.sbatch

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



