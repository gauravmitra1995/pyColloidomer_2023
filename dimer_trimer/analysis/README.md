**Analysis for dimer/trimer simulations**



**A) Information on content of folders:**

*unbondedpatches_vs_time_data* contains data for fraction of total number of free available binders as a function of simulation time for all the different simulation conditions we have explored (mainly varying number of binders (Np), droplet radii (R) and binding affinity epsilon). Each of these data files are generated from averaging over the 10 seeds run for each condition. 

*saturation_data* contains the data for the value of this fraction at saturation for various choices of R when (a) Np is varied and epsilon is fixed (b) epsilon is varied and Np is fixed. This data is basically used to construct all the 6 heat maps shown in Fig.3 of the main text and the one in Fig.S2 of SI.

*final_figures* contain all the heat maps (Fig.3, Fig.S2) and also the plots showing the variation of the fraction of available binders with time (and the fitted curves are shown too). 



**B) Instructions and sequence for running the scripts to generate the figures used in main text and SI:**

1.For each simulation condition, combine all gsd files into a single gsd file corresponding to succesive restart runs using: *bash job_combine_multipleruns*

2.Perform analysis on this combined gsd file to obtain data for the number of binders remaining on each droplet as a function of simulation time: *bash unbondedpatches_jobs_submit*

3.This will generate data for each condition separately for each of the seeds run. In order to average the data over all the seeds run for a particular simulation condition, do: *bash job_averageallseeds*

4.To generate the saturation fraction data required for plotting the heat maps as well as to do the curve fitting for the fraction vs time curves, do:

*wrapper=./../../dybond/run-hoomd2.9.6.bash*


*$wrapper python unbondedpatches_vs_time_varyRandeps.py --Nclusters 2 --Np 100* 

*$wrapper python unbondedpatches_vs_time_varyRandeps.py --Nclusters 3 --Np 100*


*$wrapper python unbondedpatches_vs_time_varyRandNp.py --Nclusters 2 --epsilon 20.7*

*$wrapper python unbondedpatches_vs_time_varyRandNp.py --Nclusters 3 --epsilon 20.7*

5.To generate the heat maps in Fig.3 and Fig.S2, do:

*wrapper=./../../dybond/run-hoomd2.9.6.bash*


*$wrapper python saturationfraction_varyingRandeps.py --Nclusters 2 --Np 100*

*$wrapper python saturationfraction_varyingRandeps.py --Nclusters 3 --Np 100*


*$wrapper python saturationfraction_varyingRandNp.py --Nclusters 2 --epsilon 20.7*

*$wrapper python saturationfraction_varyingRandNp.py --Nclusters 3 --epsilon 20.7*

6.In order to unwrap the simulation trajectory, align the droplets along a given direction (x/y/z) and bring them to the center of the simulation box, do: *bash job_align*

The final snapshots of dimers/trimer runs shown in Fig.4 of the main text could be generated after performing this step. 

**(Note: The wrapper script will work only when the singularity files are present inside the folder *./../../dybond/singularity*.)**











