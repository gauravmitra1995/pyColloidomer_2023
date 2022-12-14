#!/bin/bash

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node 1
#SBATCH --cpus-per-task 1

#SBATCH --mem 4GB

#SBATCH -t 48:00:00
#SBATCH --dependency=singleton

##SBATCH --gres=gpu:1
#SBATCH --gres=gpu:1g.10gb:1

current_dir=$(pwd)

host=$(env - hostname -s)
exe=run_simulation.py

scriptdir=$(cd $(dirname $0);pwd)
wrapper=$scriptdir/../../dybond/run-hoomd2.9.6.bash
$wrapper python -u $basedir/folding/simulation_setup/${exe}

