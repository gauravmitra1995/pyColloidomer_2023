**Analysis for folding studies**


**A) Information on content of folders:**

*folding_singletrajectory_forpaper* contains the average valence data over simulation time for the specific condition and seed we ran those 15 alternate heating/cooling cycles for. Due to file upload size constraints, the gsd file could not be kept here as well. The sub-folder *logfiles_singletrajectory* contains all the temperature and timestep information for all of the 30 runs. 

*valence_data* contains valence data for all 7 droplets over simulation time for the 300 different seeds we ran for the same condition (to quantify the populations of the folded structures obtained). Each of these data files corresponds to 5 heating/cooling cycles. 

*final_figures* contain all the figures shown in Fig.7 and Fig.S9.



**B) Instructions and sequence for running the scripts to generate the figures used in main text and SI:**

1.For each simulation condition, combine all gsd files into a single gsd file corresponding to succesive restart runs using: *bash job_combine_multipleruns*


2.To generate the temperature vs simulation time plot (also showing the variation of average valence) for the single long trajectory of 15 heating/cooling cycles (shown in Fig.7a, do:

The average valence data as well as the log files for the single long trajectory are already present in *folding_singletrajectory_forpaper*, so one can directly do as below:

*$wrapper python temperature_and_avgvalence_vs_time.py --fileprefix [....]*   

#(example fileprefix: ./folding_singletrajectory_forpaper/folding_restl2.0_Nc7_Np200_R20.0_rB1.0_kAB200.0_eps4.6_dim2_kspring10.0_gammaA0.1_gammapatch0.0001_seed12197099) 


3.To generate the histogram of populations of folded structures (Fig.7b), do:

For each of the seeds, perform the valence analysis on the combined trajectory file to get the valence data for all 7 droplets over the simulation time. 

*bash linkageanalysis_jobs_submit* 

This data will be saved in the *valence_data* folder. Then do:

*wrapper=./../../dybond/run-hoomd2.9.6.bash*

*$wrapper python histogram_structuresformed.py*


4.Finally, to generate the plot of the number of bonds of C-C and D-D types between adhesion patches of adjacent droplets with simulation time (Fig.S9),do:

*wrapper=./../../dybond/run-hoomd2.9.6.bash*

*$wrapper python bindersinpatches_analysis.py --trajectory_file [....]*

#(example trajectory_file: ./folding_singletrajectory_forpaper/folding_restl2.0_Nc7_Np200_R20.0_rB1.0_kAB200.0_eps4.6_dim2_kspring10.0_gammaA0.1_gammapatch0.0001_seed12197099.allruns.unwrap.gsd)


**(Note: The wrapper script will work only when the singularity files are present inside the folder *./../../dybond/singularity*.)**

