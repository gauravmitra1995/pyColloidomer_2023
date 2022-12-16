**Simulation setup for studying assembly of lattice of droplets with binders**


**Information on the files:**

This folder will create folders for each of the conditions we have and run MD simulations inside those folders. All relevant simulation parameters are inside the files *input_general_lattice.yaml*, *input_clusters_lattice.yaml* and *input_particles.yaml* (also mentioned in Tables S3 to S6 in the paper).

The bash script *vary_parameters_lattice.sh* and *continuesim.sh* takes in important parameters which need to be provided from outside by the user and overrides the default values in the yaml files, using the *update_yaml_lattice.py* script. 

*run-all.sbatch* is the sbatch script for submitting jobes in queue to the HPC cluster.

*finalframes* is a folder containing the final frames from the combined gsd files (all the successive runs for a given simulation condition). The final frame files ending with "coloring_by_structure.gsd" were used to prepare the snapshots shown in Fig. 5 of the main text and supplemental figures S6 and S7.


**Instructions:**

1.Copy the main simulation run scripts (5 of them) from the folder main_python_scripts_simulationsetup to here

2.To run a simulation for the first time, do *bash vary_parameters_lattice.sh*
  For any successive restart run, do *bash continuesim.sh*

3.All the simulation directories will be created and the updated yaml files and sbatch script also get copied to those directories. 


     
