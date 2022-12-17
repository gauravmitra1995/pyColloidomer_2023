This folder contains the code for dynamic binding / unbinding between particles. 

**(Singularity files have to be present inside the folder *dybond/singularity/* for the compilation of the plugin to work properly!!!!)**


**Information about important files in this folder**

The folder *nvidia-hoomd-2.9.6* contains the script *DyBondUpdater.cc* which basically contains the code for this dynamic binding. It also contains a header file *DyBondUpdater.h* where all the important global variables, functions and constructors are defined (which are used in the .cc script). A python script *update.py* is also present which activates the c++ DyBondUpdater from cpp module by initializing the dybond plugin and also feeds parameters to the DyBondUpdater written in C++ which are sent from the python simulation run script. 


**Instructions for compilation of the dybond plugin (with singularity present)**

1.Go to the folder *nvidia-hoomd-2.9.6*.
2.Then do: *./run-hoomd-2.9.6.bash*
3.Do: *bash compile.bash*
4.If the dybond plugin compiles correctly, a file named *_dybond_plugin.so* will be generated inside *nvidia-hoomd-2.9.6*.
5.The wrapper script to finally use for running simulations or performing analyses is: *run-hoomd2.9.6.bash* present inside the folder *dybond*.



