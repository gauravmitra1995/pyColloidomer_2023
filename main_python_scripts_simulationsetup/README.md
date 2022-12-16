**Python framework for the CG model: setup of the system and MD simulations**

This folder contains 5 scripts:

1.*read_parameters.py* reads in  parameters from the yaml files, defines important functions to assign important system parameters.

2.*interactions.py* defines a function to setup the tabulated soft quartic potential in HOOMD.

3.*clustermaker.py* designs each droplet coated with droplets, assigns the positions and physical attributes of the binders, creates bonds and angle terms between droplet and binders.

4.*systemsetup.py* sets up a snapshot in HOOMD from all the system information (particle attributes, bonds/angles), implements the bonded and non-bonded interactions and also defines function to call the DynBondUpdater written as a plugin to HOOMD-Blue in C++ (scripts available in the *dybond* folder of this repository).

5.*run_simulation.py* basically invokes everything that the other 4 scripts establish and sets up a Langevin Dynamics simulation for this system of droplets with binders. The simulation trajectory is dumped into a .gsd file in the run directory. 
