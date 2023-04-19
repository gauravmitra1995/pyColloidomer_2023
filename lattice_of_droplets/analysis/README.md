**Analysis for simulations to study assembly of lattice of droplets with binders**




**A) Information on content of folders:**

*valence_data* constains droplet valence data (for all 81 droplets) with time (combined over all the seeds) as well as the average droplet valence data with time

*structuralanalysis_data* contains fraction of droplets in different structures with time (combined over all the seeds)

*plot_data* contains fraction of droplets in different valences with time (combined over all the seeds)

*histograms_finalframe* contains the fraction of droplets in different valences / structures in the final frame (combined over all seeds)

*finalframefractions_valences_averageallseeds* contains the average and standard deviations of the fraction of droplets in different valences in the final frame (here, averaging is done over the 10 seeds)

*finalframefractions_structures_averageallseeds* contains the average and standard deviations of the fraction of droplets in different structures in the final frame (here, averaging is done over the 10 seeds)

*chainlengths_finalframe* contains the statistics for the average and standard deviations of the maximum chain lengths in the final frame and also the fractions of various chain sizes from a combined bootstrapped data over all seeds (required for making the histogram)

*final_figures* has all the relevant figures shown in text and SI. 




**B) Instructions and sequence for running the scripts to generate the figures used in main text and SI:**

1.For each simulation condition, combine all gsd files into a single gsd file corresponding to succesive restart runs using: *bash job_combine_multipleruns*

2.Perform analysis of bond valence of the droplets and clustering of the structures obtained using: *bash analysis_jobs_submit*

3.Unwrap the trajectories for every binder particle in the system using: *bash job_unwrapping*

4.Get the bond valence data, structural clustering data as well as the histograms of bond valence and structures formed and the chain length information averaged/combined over all seeds using: 

*bash job_averageallseeds*

5.Finally obtain the bond valence vs time plots shown in Fig.6a and Fig.S10, the histogram of bond valences at final frame shown in Fig.6b and Fig.S6, the histogram of final structures shown in Fig.6c and Fig.S5, and the histogram of chain lengths shown in Fig.6d and Fig.S7, using:

*bash job_plotting*

6.To get the bar charts showing the 3 metrics for optimizing the yield of colloidomers, for the 8 different conditions of (phi,gamma), shown in Fig.5, do:

*wrapper=./../../dybond/run-hoomd2.9.6.bash*

*wrapper python plot_bars_optimization_metrics.py --Np 100 --R 50.0 --epsilon 20.7 --Nclusters 81 --radiusB 1.0 --gammapatch 0.0001 --kspring 10.0 --restlength 2.0*        

7.To get the final frames of the simulation trajectories and then unwrap them, do:

*bash job_getfinalframegsd_unwrapped*

8.Finally, to color the droplets (and their associated binders) as per the structure to which they belong / droplet valence, do:

*bash job_coloring*  

**(Note: The wrapper script will work only when the singularity files are present inside the folder *./../../dybond/singularity*.)**
