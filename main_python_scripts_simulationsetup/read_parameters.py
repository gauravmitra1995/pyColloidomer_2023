import numpy as np

def softrepulsion(r,rmin,rmax,epsilon,ron,rcut):

    if(r<ron):
        S=1
        dS=0
    elif(ron<=r and r<=rcut):
        S=((rcut**2-r**2)**2)*(rcut**2 + 2*r**2 - 3*ron**2)/((rcut**2-ron**2)**3)
        dS=((12*(r**5))-(12*(rcut**2)*(r**3))-(12*(ron**2)*(r**3))+(12*(rcut**2)*(ron**2)*r))/((rcut**2-ron**2)**3)
    elif(r>rcut):
        S=0
        dS=0

    if(r<rcut):
        U=epsilon*(1-(r/rcut)**4)
        dU=epsilon*(-4*(r**3)/(rcut**4))
    else:
        U=0
        dU=0

    shift=0

    if(ron<rcut):
        V=S*U
        F=-(S*dU + U*dS)
    elif(ron>=rcut):
        V=U-shift
        F=-dU

    return(V,F)

def spherical2cart(rho,theta,phi):
    x=rho*np.cos(theta)*np.sin(phi)
    y=rho*np.sin(theta)*np.sin(phi)
    z=rho*np.cos(phi)
    return [x,y,z]

def sphere_fibonacci_grid_points ( ng,r0AB,r0AC ):   #THE MORE AUTHENTIC ALGORITHM!!!!

  phi = ( 1.0 + np.sqrt ( 5.0 ) ) / 2.0

  theta = np.zeros ( ng )
  sphi = np.zeros ( ng )
  cphi = np.zeros ( ng )

  for i in range ( 0, ng ):
    i2 = 2 * i - ( ng - 1 )
    theta[i] = 2.0 * np.pi * float ( i2 ) / phi
    sphi[i] = float ( i2 ) / float ( ng )
    cphi[i] = np.sqrt ( float ( ng + i2 ) * float ( ng - i2 ) ) / float ( ng )

  posB = np.zeros ( ( ng, 3 ) )
  posC = np.zeros ( ( ng, 3 ) )

  for i in range ( 0, ng ) :
    posB[i,0] = cphi[i] * np.sin ( theta[i] ) * r0AB
    posB[i,1] = cphi[i] * np.cos ( theta[i] ) * r0AB
    posB[i,2] = sphi[i] * r0AB

    posC[i,0] = cphi[i] * np.sin ( theta[i] ) * r0AC
    posC[i,1] = cphi[i] * np.cos ( theta[i] ) * r0AC
    posC[i,2] = sphi[i] * r0AC

  return posB,posC

def rotation(angle,axis):
    c = np.cos(angle)
    s = np.sin(angle)
    if axis=='x':
        return np.array(
            [[1, 0, 0],
             [0, c, -s],
             [0, s, c]])
    elif axis=='y':
        return np.array(
            [[c, 0, s],
             [0, 1, 0],
             [-s, 0, c]])
    elif axis=='z':
        return np.array(
            [[c, -s, 0],
             [s, c, 0],
             [0, 0, 1]])

def calculate_distance(r1,r2):
    r=r2-r1
    dist=np.sqrt(r[0]**2+r[1]**2+r[2]**2)
    return dist

def flat_list(alist):
    flatlist=list()
    for sublist in alist:
        for i in sublist:
            flatlist.append(i)
    return(flatlist)


def load_yaml():

    ####################################################### Loading the yaml input files #################################################
    import yaml

    with open ('input_general.yaml') as f1:
        data_general = yaml.load(f1,Loader=yaml.FullLoader)
    with open ('input_clusters.yaml') as f2:
        data_clusters = yaml.load(f2,Loader=yaml.FullLoader)
    with open ('input_particles.yaml') as f3:
            data_particles = yaml.load(f3,Loader=yaml.FullLoader)

    return data_general,data_clusters,data_particles


def read_yamlfiles(data_general,data_clusters,data_particles):

    ####################################################### Reading general info ########################################################

    generalinputdict=data_general # create a dictionary for general inputs

    ####################################################### Reading clusters info #######################################################

    clustertypeinfolist= list() # a list of dictionaries, each dictionary contains info for each cluster type

    for i in data_clusters.keys():
        clusterinfodict={}
        key=i
        clusterinfodict['Np_per_sup']=data_clusters[key]['Np_per_sup']
        clusterinfodict['typelist']=key
        clusterinfodict['Nc']=data_clusters[key]['Nc']
        clustertypeinfolist.append(clusterinfodict)

    ################################################ Reading particle size and interactions ##############################################

    input_r_dict={}
    input_mass_dict={}

    softV_list=list()
    input_softV_epsilon_dict={}
    input_softV_rcut_dict={}

    input_k_dict={}
    input_r0_dict={}

    input_dybond_dict={}

    for i in data_particles['particle'].keys():
        dic=data_particles['particle'][i]
        key=i
        input_r_dict[key]=dic['radius']
        input_mass_dict[key]=dic['mass']

    for i in data_particles['harmonic'].keys():
        dic=data_particles['harmonic'][i]
        key=i
        input_k_dict[key]=dic['kspring']
        input_r0_dict[key]=dic['r0']

    for i in data_particles['soft_V'].keys():
        dic=data_particles['soft_V'][i]
        key=i
        key="".join(sorted(key))
        softV_list.append(key)
        input_softV_epsilon_dict[key]=dic['epsilon']
        input_softV_rcut_dict[key]=dic['rcut']

    for i in data_particles['dybond'].keys():
        dic=data_particles['dybond'][i]
        key=i
        input_dybond_dict[key]=dic

    return generalinputdict,clustertypeinfolist,input_r_dict,input_mass_dict,input_k_dict,input_r0_dict,softV_list,input_softV_epsilon_dict,input_softV_rcut_dict,input_dybond_dict


################################## Set up the simulation to be restartable if progressfile is found ##################################

def restart_setup(fileprefix,nsteps):
    
    import os
    progressfile = fileprefix+'.progress.txt'
    if os.path.exists(progressfile):
        continuesim=True
        prevsteps = int(open(progressfile,'r').readlines()[0])
        inputfile = fileprefix+'.run.%i.gsd'%prevsteps
        totalsteps = prevsteps + nsteps
    else:
        continuesim = False
        prevsteps = 0
        totalsteps = prevsteps + nsteps
        inputfile=None

    outputgsdfile = fileprefix+'.run.%i.gsd'%totalsteps
    outputlogfile = fileprefix+'.run.%i.log'%totalsteps

    return continuesim,totalsteps,outputgsdfile,outputlogfile,progressfile,inputfile


########################################################### Getting list of all particle types ############################################################
########################################################### Getting a dictionary of particle types for each cluster and their corresponding Np's ############################################################


def get_info_droplet(clustertypeinfolist):
    
    type_list = list()
    terminal_particles = list()

    for i in range(len(clustertypeinfolist)): # first index in loop: cluster type
        tlist=[t for t in clustertypeinfolist[i]['typelist']] # list of types for each clustertype
        terminal_particles.extend(tlist[2:len(tlist)])  #Since 1st two types in the type list are A and B which are same in every cluster

    type_list=tlist[:2] + terminal_particles
    terminal_particles=sorted(set(terminal_particles))

    cluster_particletypedict={}
    Np_dict={}
    Nparticles=0
    clusid=0

    for i in range(len(clustertypeinfolist)): # first index in loop: cluster type
        tlist=[t for t in clustertypeinfolist[i]['typelist']] # list of types for each clustertype
        clusterids=[]
        if(type(clustertypeinfolist[i]['Nc']) is list):
            clusterids=clustertypeinfolist[i]['Nc']
            for j in range(len(clusterids)): # second index in loop: cluster
                if(type(clustertypeinfolist[i]['Np_per_sup']) is list):
                    Np_per_sup=clustertypeinfolist[i]['Np_per_sup'][j]
                else:
                    Np_per_sup=clustertypeinfolist[i]['Np_per_sup']
                Np_dict[clusterids[j]]=Np_per_sup
                cluster_particletypedict[clusterids[j]]=tlist

        else:
            for j in range(clustertypeinfolist[i]['Nc']): # second index in loop: cluster
                Np_dict[clusid]=clustertypeinfolist[i]['Np_per_sup']
                cluster_particletypedict[clusid]=tlist
                clusid+=1

    for i in Np_dict.keys():
        Np=Np_dict[i]
        Nparticles+=2*Np+1

    Nclus=len(list(Np_dict.keys()))

    return type_list,Np_dict,cluster_particletypedict,Nparticles,Nclus


########################################################### setting radii and masses ############################################################


def set_radii_masses(type_list,input_r_dict,input_mass_dict):

    r_dict={}
    mass_dict={}

    for i in range(len(type_list)):
        if type_list[i] in input_r_dict.keys():
            r_dict[type_list[i]]=float(input_r_dict[type_list[i]])
        else:
            r_dict[type_list[i]]= float(min(list(input_r_dict.values())))     #default radius set to smallest input value

    for i in range(len(type_list)):
        if type_list[i] in input_mass_dict.keys():
            mass_dict[type_list[i]]=float(input_mass_dict[type_list[i]])
        else:
            mass_dict[type_list[i]]=1.0     #default mass

    return r_dict,mass_dict


############################# Calculating dynamic bond length and populating the r0 dict for harmonic bonds ##########################

def set_harmonic_r0(zerodynbondlength,type_list,r_dict,input_r0_dict):

    if(zerodynbondlength=='True'):
        r0=0.0
    else:
        r0=2.0*r_dict[type_list[2]]

    for key,value in input_r0_dict.items():
        if value is None:
            if(len(key)==2):
                if(key[0]=='A' and key[1]=='B'):
                    input_r0_dict[key]=r_dict['A']+r_dict['B']
                elif(key[0]=='B'):
                    input_r0_dict[key]=r_dict['B']+r_dict['C']
            else:
                input_r0_dict[key]=np.pi

    for key,value in input_r0_dict.items():
        if value is None and key[0]=='A' and key[1]!='B':
            input_r0_dict[key]=input_r0_dict[key[0]+'B']+input_r0_dict['B'+key[1]]

    return r0,input_r0_dict
                                                                


############################### Specifying the lattice constant and box_dimensions for the respective type of simulation ###############################


def get_boxsize(chainlink,areafraction,lattice_dim,Nclus,r_dict,type_list,dimension):
    
    if(chainlink=='False' and areafraction is not None):
        dim1=lattice_dim[0]
        dim2=dim1
        squarelc=(Nclus*3.142*(r_dict[type_list[0]]**2))/(areafraction*(dim1)*(dim1))
        latticeconst=float(format(np.sqrt(squarelc),'.2f'))
        Lx=((lattice_dim[0])*latticeconst)     #Defining Lx,Ly,Lz for cluster arrangement into a lattice; not the actual simulation box dimension
        Ly=((lattice_dim[1])*latticeconst)
        Lz=dim1*2.0*r_dict[type_list[0]]
    else:
        dim1=lattice_dim[1]
        latticeconst=2.0*(r_dict[type_list[0]]+2.0*r_dict[type_list[1]]+r_dict[type_list[2]]) + 2.0*r_dict[type_list[2]]
        maxrcut=dim1*2.0*r_dict[type_list[0]]
        Lx=5*maxrcut
        Ly=5*maxrcut
        Lz=3*maxrcut
    
    #if linear chain (one of the box dimensions is 1), then expand the box to a cube (Lx=Ly=Lz) for 3-D and for 2-D, make Lx and Ly equal and not equal to Lz

    if lattice_dim[0]==1 and lattice_dim[1]!=1:
         Lx=Ly
         if(dimension==3):
             Lz=Ly
    elif lattice_dim[1]==1 and lattice_dim[0]!=1:
         Ly=Lx
         if(dimension==3):
             Lz=Lx
    elif lattice_dim[0]==1 and lattice_dim[1]==1:
         Lx=Ly
         if(dimension==3):
             Lz=Ly

    box_size=np.array([Lx,Ly,Lz])

    return latticeconst,box_size

########################################### Create lattice sites) for central particles to be located and arrange them (alternate or random) ########################################

def get_latticepoints(latticeconst,lattice_dim,sequence): 
   
    import random 

    Lx=((lattice_dim[0])*latticeconst)
    Ly=((lattice_dim[1])*latticeconst)

    latticepoints=list()
    x0=-Lx/2 + latticeconst/2
    y0=Ly/2 - latticeconst/2

    for i in range(lattice_dim[0]):
        for j in range(lattice_dim[1]):
            latticepoints.append([(x0+latticeconst*i),(y0-latticeconst*j),0])

    latticepoints.sort(key=lambda x: x[1])  #sorting the lattice points in the right order, according to their position

    if sequence=='alternate': #[0,1,2,3,4,5]-->[0,2,4,1,3,5]; this new list of lattice points, like the old one, is sequentially assigned to each newly created cluster object
         latticepoints=[latticepoints[i] for i in range(len(latticepoints)) if i%2==0]+[latticepoints[i] for i in range(len(latticepoints)) if i%2==1]
    elif sequence=='random': # then arrange randomly
         latticepoints=random.sample(latticepoints,len(latticepoints))

    return latticepoints


###################################################### Setting parameters for dynamic bonds #################################################

def set_dybondparameters(input_dybond_dict,r0,Nparticles,Np_dict,Nclus):

    particle=0
    clusterids_allparticles=np.zeros(Nparticles)
    for i in range(Nclus):
        Np=Np_dict[i]
        for j in range(particle,particle+(2*Np+1),1):
            clusterids_allparticles[j]=i
            particle=particle+1

    dybond_dict={}
    for k in list(input_dybond_dict.keys()):
        dybond_dict[k[0]+','+k[1]]=input_dybond_dict[k]
        if(dybond_dict[k[0]+','+k[1]]['r0'] is None):
            dybond_dict[k[0]+','+k[1]]['r0']=r0

    particle_clusterids=[]
    for i in range(Nparticles):
        particle_clusterids.append(list(clusterids_allparticles)[i])
    particle_clusterids=' '.join(str(i) for i in particle_clusterids)

    return dybond_dict,particle_clusterids

############################################# Setting soft V parameters (epsilon, rcut) ##############################################


def set_softVparameters(softV_list,input_softV_epsilon_dict,input_softV_rcut_dict,input_r0_dict,r_dict):

    softV_epsilon_dict={}
    softV_rcut_dict={}

    for i in softV_list:
        if(input_softV_epsilon_dict[i[0]+i[1]]!=None):
            softV_epsilon_dict[i[0]+','+i[1]]=float(input_softV_epsilon_dict[i[0]+i[1]])
        else:
            softV_epsilon_dict[i[0]+','+i[1]]=200.0

        if(input_softV_rcut_dict[i[0]+i[1]]!=None):
            softV_rcut_dict[i[0]+','+i[1]]=float(input_softV_rcut_dict[i[0]+i[1]])
        else:
            if(i[0]=='A' and i[1]!='B' and i[1]!='A'):
                softV_rcut_dict[i[0]+','+i[1]]=float(1.0*(input_r0_dict[i[0]+i[1]]))
            elif(i[0]=='A' and i[1]=='A'):
                softV_rcut_dict[i[0]+','+i[1]]=float(1.1*2.0*r_dict[i[0]])
            else:
                softV_rcut_dict[i[0]+','+i[1]]=float(1.0*(r_dict[i[0]]+r_dict[i[1]]))

    return softV_epsilon_dict,softV_rcut_dict


############################################# Setting wall parameters ##############################################


def set_wall_parameters(dimension,r_dict,type_list):

    up_wall_z=1.25*2.0*r_dict[type_list[0]]
    low_wall_z=-1.25*2.0*r_dict[type_list[0]]

    if(dimension==2):   #specifying wall_epsilon based on the dimensionality of system
        wall_epsilon=10.0
    else:
        wall_epsilon=0

    return wall_epsilon,up_wall_z,low_wall_z 
