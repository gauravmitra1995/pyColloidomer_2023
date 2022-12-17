**A Coarse-Grained Simulation Model for Self-Assembly of Liquid Droplets Featuring Explicit Mobile Binders (SUBMITTED)**

Authors:

**Gaurav Mitra, Chuan Chang, Angus McMullen, Daniela Puchall, Jasna Brujic, and Glen M. Hocky**


This repository contains the python framework for performing MD simulations to study self-assembly of droplets with mobile binders.

We have three main folders/directories for the 3 types of simulations we perform:

1.dimer_trimer
2.lattice_of_droplets
3.folding

The dynamic binding/unbinding code is in the folder *dybond*. For successful compilation of the dybond plugin to HOOMD-Blue as well as for the wrapper script to work,

**Singularity files have to be present in a folder *singularity* inside ./dybond/**

Each folder for the 3 types of simulations contains an *analysis* folder and a *simulation_setup* folder. In addition, there is a folder named *main_python_scripts_simulationsetup* outside of all these folders which contains the 5 main .py scripts which can read in parameters from yaml files, set up the system and run the MD simulations by calling the Dynamic Bond Updater.

Detailed instructions for setting up simulations and performing specific analyses (and generating the figures discussed in the main text / SI) is provided in the README files of these sub-folders. 


