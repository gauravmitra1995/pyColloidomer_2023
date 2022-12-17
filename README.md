#Manuscript: **A Coarse-Grained Simulation Model for Self-Assembly of Liquid Droplets Featuring Explicit Mobile Binders**  **(SUBMITTED)**

#Authors: **Gaurav Mitra, Chuan Chang, Angus McMullen, Daniela Puchall, Jasna Brujic, and Glen M. Hocky**


This repository contains the python framework for performing MD simulations to study self-assembly of droplets with mobile binders.

We have three main folders/directories for the 3 types of simulations we perform:

1.dimer_trimer
2.lattice_of_droplets
3.folding

The dynamic bonding code is in the folder "dybond".

Each folder for the 3 types of simulations contains an "analysis" folder and a "simulation_setup" folder. In addition, there is a folder named "main_python_scripts_simulationsetup" outside of all these folders which contains the 5 main .py scripts which can read in parameters from yaml files, set up the system and run the MD simulations by calling the DynBond Updater. These 5 scripts have to be first copied into the "simulation_setup" folder before setting up simulations i.e. do the following:

Inside the folder pyColloidomer_Gaurav:
cd main_python_scripts_simulationsetup
cp read_parameters.py ../dimer_trimer 
cp run_simulation.py ../dimer_trimer 
.
.
. 
and so on.......
