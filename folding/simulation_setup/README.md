**Simulation setup for folding studies for chains**

**(Assuming singularity is present inside ./../../dybond/)**


**A) Information on the files:**

A hierarchy of folders will be created for each unique combination of simulation parameters we have and MD simulations will run inside the final folder of that hierarchy. All relevant simulation parameters are inside the files *input_general.yaml*, *input_clusters.yaml* and *input_particles.yaml* (also mentioned in Tables S3,S6 in the paper).

The bash script *vary_parameters_folding.sh* and *continuesim.sh* takes in important parameters which need to be provided from outside by the user and overrides the default values in the yaml files, using the *update_yaml_folding.py* script. The bash scripts also help set up the run directories.

*run-all.sbatch* is the sbatch script for submitting jobs in queue to the HPC cluster.


**B) Instructions:**

1.Copy the main simulation run scripts (5 of them) from the folder *main_python_scripts_simulationsetup* to here.

2.Change/update the parameters inside the bash scripts (*vary_parameters_lattice.sh* or *continuesim.sh*) if needed. 

3.If simulation is to be run on CPU, set "gpu" to 0 in the *input_general.yaml* file, otherwise for gpu runs, set it to 1. Also check the values for no of steps to run: "num_step", timestep: "dt" and other such relevant parameters in the yaml files before submitting a job. 

4.Modify the path of the base directory *pyColloidomer_2023* in the *run-all.sbatch* script before starting a simulation.

5.To run a simulation for the first time, do *bash vary_parameters_folding.sh*.

6.For any successive restart run, do *bash continuesim.sh*.

All the simulation directories will be created and the updated yaml files and sbatch script also get copied to those directories. 



