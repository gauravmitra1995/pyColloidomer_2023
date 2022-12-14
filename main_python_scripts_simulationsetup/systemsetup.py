import numpy as np

############################################### Functions to create droplets with patches, assign bonds and angles ###########################################


def create_clusters(Np_dict,cluster_particletypedict,latticepoints,r_dict,mass_dict,input_r0_dict,choice_initialdistribution,chainlink):
    
    from clustermaker import patchyparticlecluster

    objectlist=list()
    totalNsofar=0
    totalNbonds=0
    totalNangles=0

    print("Typelist,  Clusterid,  Np,  Centerpos")
    for i in range(len(Np_dict.keys())):
        clusid=list(Np_dict.keys())[i]
        tlist=cluster_particletypedict[clusid]
        cluster=patchyparticlecluster(Np_per_sup=Np_dict[clusid],typelist=tlist,centerpos=latticepoints[clusid],clusterid=clusid,r_dict=r_dict,mass_dict=mass_dict,input_r0_dict=input_r0_dict,choice_initialdistribution=choice_initialdistribution)
        objectlist.append(cluster)

    objectlist=sorted(objectlist,key=lambda x:x.centerpos)
    
    for i in range(len(objectlist)):
        cluster=objectlist[i]
        cluster.clusterid=i
        print(cluster.typelist,cluster.clusterid,cluster.Np,cluster.centerpos)
        cluster.assign_patches(chainlink=chainlink)
        cluster.set_bonds_angles(previousN=totalNsofar)
        totalNsofar+=cluster.N
        totalNbonds+=cluster.Np*2
        totalNangles+=cluster.Np

    return objectlist,totalNbonds,totalNangles



def assign_particleproperties(snapshot,objectlist,type_list,box_size):
    particle=0
    for i in range(len(objectlist)):
        cluster=objectlist[i]
        typeids=[type_list.index(p) for p in cluster.particle_types]
        cluster.particle_typeids=typeids
        cluster.particle_positions=[xyz-box_size*(np.floor((xyz/box_size)+0.5)) for xyz in cluster.particle_positions]
        for x in range(len(cluster.particle_positions)):
            snapshot.particles.position[particle]=cluster.particle_positions[x]
            snapshot.particles.diameter[particle]=cluster.particle_diameters[x]
            snapshot.particles.mass[particle]=cluster.particle_masses[x]
            snapshot.particles.typeid[particle]=cluster.particle_typeids[x]
            particle+=1



def combine_bondsangles_allparticles(objectlist):

    all_harmonic_bonds_dict={} # key: bond pair (e.g. AB), value: list of bond pairs, collected from all cluster objects that have this bond pair
    all_harmonic_angles_dict={} # key: angle pair (e.g. ABC,), value: list of angle pairs, collected from all cluster objects that have this angle pair
    bondtypes=[]
    angletypes=[]

    for i in range(len(objectlist)): # loop over all clusters
        cluster=objectlist[i]

        bondtypes.extend(k for k in list(cluster.bonds.keys()) if k not in bondtypes)
        angletypes.extend(k for k in list(cluster.angles.keys()) if k not in angletypes)

        for bond_type in cluster.bonds.keys(): # loop over the bond types of cluster i
            if bond_type not in all_harmonic_bonds_dict.keys():
                all_harmonic_bonds_dict[bond_type]=cluster.bonds[bond_type]
            else:
                all_harmonic_bonds_dict[bond_type].extend(cluster.bonds[bond_type])

        for angle_type in cluster.angles.keys(): #loop over the angle types of cluster i
            if angle_type not in all_harmonic_angles_dict.keys():
                all_harmonic_angles_dict[angle_type]=cluster.angles[angle_type]
            else:
                all_harmonic_angles_dict[angle_type].extend(cluster.angles[angle_type])

    return all_harmonic_bonds_dict, all_harmonic_angles_dict,bondtypes,angletypes



def define_snapshotbondsangles(snapshot,bondtypes,angletypes,totalNbonds,totalNangles,input_dybond_dict,all_harmonic_bonds_dict,all_harmonic_angles_dict):

    snapshot.bonds.resize(totalNbonds)
    snapshot.angles.resize(totalNangles)
    snapshot.bonds.types=bondtypes+[(pair[0]+'-'+pair[1]) for pair in input_dybond_dict.keys()]
    snapshot.angles.types=angletypes

    bondscreated=0
    anglescreated=0

    for bond_type in all_harmonic_bonds_dict.keys():
        if len(all_harmonic_bonds_dict[bond_type])!=0:
            snapshot.bonds.group[bondscreated:bondscreated+len(all_harmonic_bonds_dict[bond_type])] = all_harmonic_bonds_dict[bond_type]
            snapshot.bonds.typeid[bondscreated:bondscreated+len(all_harmonic_bonds_dict[bond_type])] = snapshot.bonds.types.index(bond_type)
            bondscreated+=len(all_harmonic_bonds_dict[bond_type])

    for angle_type in all_harmonic_angles_dict.keys():
        if len(all_harmonic_angles_dict[angle_type])!=0:
            snapshot.angles.group[anglescreated:anglescreated+len(all_harmonic_angles_dict[angle_type])] = all_harmonic_angles_dict[angle_type]
            snapshot.angles.typeid[anglescreated:anglescreated+len(all_harmonic_angles_dict[angle_type])] = snapshot.angles.types.index(angle_type)
            anglescreated+=len(all_harmonic_angles_dict[angle_type])


def get_chainlink_ids(snapshot,type_list):
    
    terminal_particles=type_list[2:]
    indices=[i for i in range(snapshot.particles.N) if type_list[snapshot.particles.typeid[i]] in terminal_particles]
    chainlinkids=[]
    for i in range(0,len(indices)):
        if((indices[i]-indices[i-1])!=1):
            chainlinkids.append(indices[i])
            chainlinkids.append(indices[i+1])

    chainlink_ids_array=np.array(chainlinkids)
    return chainlink_ids_array


###################################################### Set harmonic bond and angle coeffs ######################################################

def set_bond_coeffs(harmonicbond,snapshot,input_k_dict,input_r0_dict,dybond_dict):

    for pair in list(snapshot.bonds.types):
        if(pair not in input_k_dict.keys()):
            splits=pair.split('-')
            key=splits[0]+','+splits[1]
            harmonicbond.bond_coeff.set(pair, k=dybond_dict[key]['kspring'], r0=dybond_dict[key]['r0'])
        else:
            harmonicbond.bond_coeff.set(pair, k=input_k_dict[pair], r0=input_r0_dict[pair])


def set_angle_coeffs(harmonicangle,snapshot,input_k_dict,input_r0_dict):

    for pair in list(snapshot.angles.types):
        harmonicangle.angle_coeff.set(pair, k=input_k_dict[pair], t0=input_r0_dict[pair])


################################################ Set wall coefficients #################################################

def set_wall_coeffs(wall_epsilon,r_dict,type_list,upper_wall_force,lower_wall_force):

    for type in type_list:
        if type==type_list[0]: # wall only affects central particles
            upper_wall_force.force_coeff.set(type,epsilon=wall_epsilon,sigma=2.0*r_dict[type_list[0]],alpha=1.0)
            lower_wall_force.force_coeff.set(type,epsilon=wall_epsilon,sigma=2.0*r_dict[type_list[0]],alpha=1.0)
        else:
            upper_wall_force.force_coeff.set(type,epsilon=0,sigma=0,alpha=1.0,r_cut=-1)
            lower_wall_force.force_coeff.set(type,epsilon=0,sigma=0,alpha=1.0,r_cut=-1)


################################################ Call the Dynamic Bond Updater ################################################ 

def call_dybond_updater(objectlist,nl,groupall,dybond_dict,particle_clusterids,totalNbonds,dt,seed):

    #import dybond plugin & run only if dybond_dict is not empty and there is more than one droplet and there are least one pair of patches in each
    
    for obj in objectlist:
        if obj.Np_per_sup >=1:
            flag=1
        else:
            flag=0
            break

    if bool(dybond_dict) and len(objectlist)>1 and flag==1:

            import hoomd.dybond_plugin as db

            updaterdict={}
            kT_init=1.0

            for key in list(dybond_dict.keys()): # 'C,D'
                pair=key.strip().split(',') # 'C,D'-->['C','D']
                bond_type='-'.join([i for i in pair]) # ['C','D']-->'C-D'
                standarddev=(np.sqrt(kT_init/dybond_dict[key]['kspring']))
                r0=float(dybond_dict[key]['r0'])

                # rmin/rmax either user-given or calculated:

                if dybond_dict[key]['rmin'] is None:
                    if(r0==0.0):
                        rmin=0.0
                    else:
                        rmin=r0-(2.0*standarddev)
                else:
                    rmin=float(dybond_dict[key]['rmin'])

                if dybond_dict[key]['rmax']is None:
                    rmax=r0+(2.0*standarddev)
                else:
                    rmax=float(dybond_dict[key]['rmax'])

                dybondchecksteps=int(dybond_dict[key]['dybondchecksteps'])
                self_avoiding=dybond_dict[key]['self_avoid_chain']
                kspring=float(dybond_dict[key]['kspring'])
                flag_temperature=dybond_dict[key]['temp_dep']
                flag_force=dybond_dict[key]['force_dep']
                flag_MP=int(dybond_dict[key]['metropolis'])
                Tmelt=float(dybond_dict[key]['Tmelt'])
                alpha=float(dybond_dict[key]['alpha'])
                kon_init=float(dybond_dict[key]['kon_init'])
                kon_melt=float(dybond_dict[key]['kon_melt'])
                koff_init=float(dybond_dict[key]['koff_init'])
                koff_melt=float(kon_init+kon_melt-koff_init)
                
                updaterdict[key] = db.update.dybond(nl,group=groupall,period=dybondchecksteps) # create a new updater for every dybond pair
                print("------------------------------------------------------------------------------")
                print("Dybond type: ",bond_type)
                print("------------------------------------------------------------------------------")
                print("Dynbond kspring: ",kspring)
                print("Dynbond rest length: ",r0)
                print("rmin: ",rmin)
                print("rmax: ",rmax)
                print("Metropolis flag: ",flag_MP)
                print("Melting temperature: ",Tmelt)
                print("Self avoiding flag: ",self_avoiding)
                print("Initial unbinding rate constant: ",koff_init)
                print("Melting unbinding rate constant: ",koff_melt)
                print("Temperature dependence flag: ",flag_temperature)
                print("//////////////////////////////////////////////////////////////////////////////")

                #Temperature dependent dynamic bond
                updaterdict[key].set_params(bond_type=bond_type,A=pair[0],B=pair[1],nondybonds=totalNbonds,r0=r0,rmin=rmin,rmax=rmax,kspring=kspring,flag_temperature=flag_temperature,flag_force=flag_force,
                                            flag_metropolis=flag_MP,Tmelt=Tmelt,alpha=alpha,kon_init=kon_init,kon_melt=kon_melt,koff_init=koff_init,koff_melt=koff_melt,checksteps=dybondchecksteps,dt=dt,
                                            particle_clusterids=particle_clusterids,self_avoiding=self_avoiding,userseed=seed)


def main():
    pass

if __name__ == "__main__":
    main()

