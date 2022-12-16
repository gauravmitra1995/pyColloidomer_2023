#!/bin/bash

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash

#Very first run set up using this script

simulationtype='lattice'
ntasks=10
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

				  seedlist=""
				  for ((x=1; x<=$ntasks; x++))
		                  do
					  no=`od -vAn -N3 -tu < /dev/urandom`
					  seedlist=`echo "$seedlist ${no}"`
			          done

				  dir=$current_dir/$simulationtype/gammaA${gammaA}_gammapatch${gammapatch}/Nclusters${Nclusters}/Np${Np}/R${R}/radiusB${rB}/areafraction${areafraction}/epsilon${epsilon}/kspring${kspring}

			          mkdir -p $dir
			
				  
				  if [ -e $dir/seedlist.txt ];then
					  echo -n " " $seedlist >> $dir/seedlist.txt
			          else
					  echo -n $seedlist > $dir/seedlist.txt
			          fi
                                  
                                  for seed in $seedlist;do

					        for zerodynbondlength in 'False';do

							rundir=$($wrapper python -u update_yaml_lattice.py --Np $Np --R $R --radiusB $rB --simulationtype $simulationtype --Nclusters $Nclusters --areafraction $areafraction --seed $seed --koninit $koninit --koffinit $koffinit --gammaA $gammaA --gammapatch $gammapatch --kspring $kspring --metropolis $metropolis --zerodynbondlength $zerodynbondlength | tail -n 1)

							echo "rundir=$rundir"

							jobname=${simulationtype}runs_gammaA${gammaA}_gammapatch${gammapatch}_Nclus${Nclusters}_Np${Np}_R${R}_rB${rB}_areafrac${areafraction}_eps${epsilon}_seed${seed}
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

