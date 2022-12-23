import numpy as np
import sys

##################################################################### Read in all the .yaml files #################################################################

from read_parameters import load_yaml,read_yamlfiles

data_general,data_clusters,data_particles = load_yaml()
generalinputdict,clustertypeinfolist,input_r_dict,input_mass_dict,input_k_dict,input_r0_dict,softV_list,input_softV_epsilon_dict,input_softV_rcut_dict,input_dybond_dict=read_yamlfiles(data_general,data_clusters,data_particles)

########################################################### global variables (read from input_general.yaml file) ##########################################################


gpu=int(generalinputdict['gpu'])
nsteps=int(generalinputdict['num_step'])
dumptime=int(generalinputdict['dumptime'])
dt=float(generalinputdict['dt'])
dumpperiod=int(dumptime/dt)
analyzeperiod=int(generalinputdict['analyzeperiod'])
fileprefix=generalinputdict['fileprefix']

seed=int(generalinputdict['seed'])
gammaA=float(generalinputdict['gammaA'])
gammapatch=float(generalinputdict['gammapatch'])
if(generalinputdict['kT']!='variant'):
    kT=float(generalinputdict['kT'])
else:
    kT=generalinputdict['kT']

lattice_dim=list(generalinputdict['lattice_dim'])
if(generalinputdict['areafraction']!=None):
    areafraction=float(generalinputdict['areafraction'])
else:
    areafraction=None
chainlink=str(generalinputdict['chain_link'])
sequence=generalinputdict['sequence']
simulationtype=generalinputdict['simulationtype']
choice_initialdistribution=int(generalinputdict['choice_initialdistribution'])
r_buff=float(generalinputdict['r_buff'])
dimension=int(generalinputdict['dimension'])
zerodynbondlength=str(generalinputdict['zerodynbondlength'])


########################################################## Assign GPU/CPU Mode #######################################################


import hoomd 

if gpu==1:
    print('GPU mode')
    hoomd.context.initialize("--mode=gpu")
else:
    print('non-GPU mode')
    hoomd.context.initialize("--mode=cpu")


############################################### Read in all relevant parameters, lists, dictionaries etc needed to create the system ###############################################


from read_parameters import restart_setup,get_info_droplet,set_radii_masses,set_harmonic_r0

continuesim,totalsteps,outputgsdfile,outputlogfile,progressfile,inputfile = restart_setup(fileprefix,nsteps)
type_list,Np_dict,cluster_particletypedict,Nparticles,Nclus = get_info_droplet(clustertypeinfolist)
r_dict,mass_dict = set_radii_masses(type_list,input_r_dict,input_mass_dict)
r0,input_r0_dict = set_harmonic_r0(zerodynbondlength,type_list,r_dict,input_r0_dict)


from read_parameters import get_boxsize,get_latticepoints

latticeconst,box_size = get_boxsize(chainlink,areafraction,lattice_dim,Nclus,r_dict,type_list,dimension)
latticepoints = get_latticepoints(latticeconst,lattice_dim,sequence)

from read_parameters import set_dybondparameters,set_softVparameters,set_wall_parameters

dybond_dict,particle_clusterids = set_dybondparameters(input_dybond_dict,r0,Nparticles,Np_dict,Nclus)
softV_epsilon_dict,softV_rcut_dict=set_softVparameters(softV_list,input_softV_epsilon_dict,input_softV_rcut_dict,input_r0_dict,r_dict)
wall_epsilon,up_wall_z,low_wall_z=set_wall_parameters(dimension,r_dict,type_list)


############################################### Create the system of droplets with patches, implement harmonic bonds/angles and interactions ###############################################


from systemsetup import create_clusters,assign_particleproperties,combine_bondsangles_allparticles,define_snapshotbondsangles

objectlist,totalNbonds,totalNangles = create_clusters(Np_dict,cluster_particletypedict,latticepoints,r_dict,mass_dict,input_r0_dict,choice_initialdistribution,chainlink)

if(continuesim==True):
    system=hoomd.init.read_gsd(inputfile,frame=-1)
else:
    snapshot= hoomd.data.make_snapshot(N=Nparticles,box=hoomd.data.boxdim(Lx=box_size[0]*1.0, Ly=box_size[1]*1.0, Lz=box_size[2]*1.0), particle_types=type_list)
    assign_particleproperties(snapshot,objectlist,type_list,box_size)
    all_harmonic_bonds_dict,all_harmonic_angles_dict,bondtypes,angletypes= combine_bondsangles_allparticles(objectlist)
    define_snapshotbondsangles(snapshot,bondtypes,angletypes,totalNbonds,totalNangles,input_dybond_dict,all_harmonic_bonds_dict,all_harmonic_angles_dict)

    system=hoomd.init.read_snapshot(snapshot)
snapshot = system.take_snapshot(bonds=True)

######################################### Set bond and angle coefficients,implement the soft V interactions and walls, call the Dynamic Bond Updater ########################################


from systemsetup import set_bond_coeffs,set_angle_coeffs,set_wall_coeffs,call_dybond_updater
import hoomd.md

harmonicbond = hoomd.md.bond.harmonic()
harmonicangle = hoomd.md.angle.harmonic()

if(totalNbonds!=0):
    set_bond_coeffs(harmonicbond,snapshot,input_k_dict,input_r0_dict,dybond_dict)

if(totalNangles!=0):
    set_angle_coeffs(harmonicangle,snapshot,input_k_dict,input_r0_dict)


from systemsetup import get_chainlink_ids

if continuesim==False:
    if chainlink=='True':
        chainlink_ids_array=get_chainlink_ids(snapshot,type_list)
        if(chainlink=='True' and bool(dybond_dict) and chainlink_ids_array.size!=0):
            for i in range(1,chainlink_ids_array.size-1):
                if(chainlink_ids_array[i+1]-chainlink_ids_array[i]>1):
                    print("Bond between ",chainlink_ids_array[i],chainlink_ids_array[i+1])
                    if(snapshot.particles.typeid[chainlink_ids_array[i]] <= snapshot.particles.typeid[chainlink_ids_array[i+1]]):
                        bond_type=type_list[snapshot.particles.typeid[chainlink_ids_array[i]]]+'-'+type_list[snapshot.particles.typeid[chainlink_ids_array[i+1]]]
                    else:
                        bond_type=type_list[snapshot.particles.typeid[chainlink_ids_array[i+1]]]+'-'+type_list[snapshot.particles.typeid[chainlink_ids_array[i]]]
                    system.bonds.add(bond_type,chainlink_ids_array[i],chainlink_ids_array[i+1])

r_buff_true=r_buff*r_dict[type_list[2]]
nl= hoomd.md.nlist.tree(r_buff=r_buff_true,check_period=1)
nl.reset_exclusions([])

table = hoomd.md.pair.table(width=1000, nlist=nl)

from read_parameters import softrepulsion
from interactions import define_interactions 

################################################ Implement soft potential #################################################


define_interactions(table,nl,softrepulsion,type_list,softV_epsilon_dict,softV_rcut_dict,r_dict)


if(wall_epsilon!=0):
    upper_wall = hoomd.md.wall.plane(origin=(0,0,up_wall_z),normal=(0,0,-1),inside=True)
    lower_wall = hoomd.md.wall.plane(origin=(0,0,low_wall_z),normal=(0,0,1),inside=True)
    wall_group = hoomd.md.wall.group(upper_wall,lower_wall)
    upper_wall_force = hoomd.md.wall.force_shifted_lj(wall_group,r_cut=(2.0**(1./6))*2.0*r_dict[type_list[0]])
    lower_wall_force = hoomd.md.wall.force_shifted_lj(wall_group,r_cut=(2.0**(1./6))*2.0*r_dict[type_list[0]])
    set_wall_coeffs(wall_epsilon,r_dict,type_list,upper_wall_force,lower_wall_force)


groupall=hoomd.group.all()
call_dybond_updater(objectlist,nl,groupall,dybond_dict,particle_clusterids,totalNbonds,dt,seed)

############################################### Run the Langevin Dynamics simulation and write to the gsd file ###############################################

hoomd.md.integrate.mode_standard(dt=dt)
ld = hoomd.md.integrate.langevin(group=groupall, kT=kT, seed=seed)

for i in type_list:
    if i=='A':
        ld.set_gamma(i,gammaA)
    else:
        ld.set_gamma(i,gammapatch)

hoomd.dump.gsd(filename=outputgsdfile, period=dumpperiod, group=groupall, overwrite=True, dynamic=['attribute','momentum','topology'])
logger=hoomd.analyze.log(filename=outputlogfile,quantities=['potential_energy', 'temperature'],period=analyzeperiod,overwrite=True)

hoomd.run(nsteps)

print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
print("No of bonds at the end: ",len(system.bonds))
bondtypeid_list=[]
for i in range(len(system.bonds)):
    typeid=system.bonds[i].typeid
    btype=system.bonds[i].type
    bondtypeid_list.append(typeid)

print(bondtypeid_list)

#very end of simulation, assuming successful end
fh = open(progressfile,'w')
fh.write("%i\n"%totalsteps)
fh.close()

