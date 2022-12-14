from functions_cluster_analysis import *
from mpl_toolkits.mplot3d import Axes3D
import random
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count


def get_clustering(cluster_bondtable_new):
    topo=cluster_bondtable_new
    l=len(topo)
    q=[] 

    q=topo

    output = []
    while len(q)>0:
        first, *rest = q
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        output.append(sorted(list(first)))
        q = rest
    
    return(output)


parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
parser.add_argument("--minframe",default=0,type=int)
parser.add_argument("--dt",default=0.001,type=float)
parser.add_argument("--frameinterval",default=1,type=int)
args = parser.parse_args()
locals().update(vars(args))

trajectory = gsd.open(trajectory_file,'rb') # read gsd file
print('file:',trajectory_file)
print('\n')

num_frames=len(trajectory)
print('Total number of frames')
print(num_frames)

maxframe=num_frames-1
######################################## Finding last uncorrupted/readable frame to use as maxframe #############################################

for i in reversed(range(num_frames)):
    print('trying frame',i)
    try:
        last_frame_testpos=trajectory[i].particles.position[0]
        last_frame=i
        break
    except:
        print('next frame')

if last_frame < maxframe:
    print('the given maxframe is unreadable, using the last readable frame, frame',last_frame,'as maxframe')
    maxframe=last_frame

##################################################################################################################################################
print('minframe used:',minframe)
print('maxframe used:',maxframe)
print('frameinterval used:',frameinterval)
frames_list=np.arange(minframe,maxframe+1,frameinterval) # list of frames to analyze
print("No of frames used for analysis:",len(frames_list))

####################################################### Gathering basic info #######################################################

system = trajectory[0] # gather basic information from 0th frame

box_size = list(system.configuration.box[:3])

type_list = list(system.particles.types) # list of all particle types in the system
num_particles = system.particles.N # total num of particles in the system


A_tags_list=[i for i in range(num_particles) if system.particles.typeid[i]==0] # tags of central particles
num_clusters=len(A_tags_list) # number of droplets
cluster_ids = np.arange(num_clusters,dtype=int)

typeid_list = [list(system.particles.typeid[i : j]) for i, j in zip([None]+A_tags_list, A_tags_list+[None]) if i!=None] # [[..cluster0 typeids..],[..cluster1 typeids..],..]

Np_list=[(len(i)-1)/2 for i in typeid_list] # list of Np's for each cluster
terminal_typeids=[i[-1] for i in typeid_list] # type ids of terminal particles (based on last typeid in each typeid sublist in typeid_list)

terminal_types=[type_list[i] for i in terminal_typeids]
terminal_types_unindexed=[i[0] for i in terminal_types]

total_clusterharmonicbonds=int(sum([i*2 for i in Np_list]))

num_non_dybonds=total_clusterharmonicbonds

bond_types=list(system.bonds.types)

dybond_types=[i for i in system.bonds.types if ((i[0] in terminal_types_unindexed) and (i[1] in terminal_types_unindexed)) or ('-' in i) or (',') in i] # if the bond type name has '-' in it, it's a dybond (maybe use some other criteria)

dybond_typeids=[system.bonds.types.index(i) for i in dybond_types]


print('Registering bonded cluster pairs and performing cluster analysis of various structures for each frame...')

structures=['monomers','dimers','linear chains','loops','other']

droplets_structuretype={}

for frame in frames_list:

    system=trajectory[int(frame)]
    bondtable=system.bonds.group.astype(int)
    cluster_bondtable=np.zeros(bondtable.shape)
    for i in range(bondtable.shape[0]):
        bond=bondtable[i]
        particle1=bond[0]
        particle2=bond[1]
        cluster1=find_cluster1(A_tags_list,particle1)
        cluster2=find_cluster1(A_tags_list,particle2)
        cluster_bondtable[i]=[int(cluster1),int(cluster2)]

    cluster_bondtable,countrepeats=remove_repeats(cluster_bondtable.tolist())
    cluster_bondtable_new=[]
    for ele in cluster_bondtable:
    	if(ele[0]!=ele[1]):
       		cluster_bondtable_new.append(ele)

    clustering=get_clustering(cluster_bondtable_new)

    droplets_in_clusters=sorted(flat_list(clustering))
    
    monomers=[]
    for x in cluster_ids:
        if(x not in droplets_in_clusters):
            monomers.append(x)

    clustering_with_singletonsadded=[]
    for i in range(len(clustering)):
        clustering_with_singletonsadded.append(clustering[i])
    for i in range(len(monomers)):
        singleton=[]
        singleton.append(monomers[i])
        clustering_with_singletonsadded.append(singleton)

    number_bonds=np.zeros(num_clusters)

    #finding the valence for every droplets and storing them in an array 
    for no in range(len(cluster_ids)):
        id=cluster_ids[no]
        for x in range(len(cluster_bondtable_new)):
            clusterno1=cluster_bondtable_new[x][0]
            clusterno2=cluster_bondtable_new[x][1]
            if(id==clusterno1 or id==clusterno2):
                 number_bonds[id]=number_bonds[id]+1

    number_of_droplets_as_monomers=0
    number_of_droplets_in_linear_chains_not_dimers=0
    number_of_droplets_in_dimers=0
    number_of_droplets_in_loops=0
    number_of_droplets_in_gelsorbranchedchains=0
    
    #assignment of droplets to the type of structure they belong to 
    for cluster in clustering_with_singletonsadded:

        valences=[]
        for i in cluster:
            valence_i=int(number_bonds[int(i)])
            valences.append(valence_i)

        no_of_1s=countX(valences,1)
        no_of_2s=countX(valences,2)

        if(len(cluster)==1):
            number_of_droplets_as_monomers+=len(cluster)
            structuretype=0
        elif(len(cluster)==2):
            if(no_of_1s==len(cluster)):
                number_of_droplets_in_dimers+=len(cluster)
                structuretype=1
        elif(len(cluster)>2):
            if(no_of_1s==0 and no_of_2s==len(cluster)):
                number_of_droplets_in_loops+=len(cluster)
                structuretype=3
            elif(no_of_1s==2 and no_of_2s==len(cluster)-2):
                number_of_droplets_in_linear_chains_not_dimers+=len(cluster)
                structuretype=2
            else:
                number_of_droplets_in_gelsorbranchedchains+=len(cluster)
                structuretype=4

        for i in cluster:
            droplets_structuretype[int(i)]=structuretype

droplets_structuretype=dict(sorted(droplets_structuretype.items(), key=lambda x:x[0]))

droplets_allparticles=[]
for particle in range(system.particles.N):
    droplet=find_cluster1(A_tags_list,particle)
    droplets_allparticles.append(droplet)

#Create new particle types based on the kind of structures to which droplets belong
system.particles.types=['P','Q','R','S','T']

particles_structuretype={}
structuretypes=np.zeros(system.particles.N)

#Label particles as per by their structure (binders assigned to a droplet get the same particle type as their droplets do)
for j in range(len(droplets_allparticles)):
    droplet=droplets_allparticles[j]
    particle=j
    for key in droplets_structuretype.keys():
        if(key==droplet):
            particles_structuretype[particle]=droplets_structuretype[droplet]
            structuretypes[particle]=int(droplets_structuretype[droplet])
            system.particles.typeid[j]=structuretypes[particle]

basename=os.path.basename(trajectory_file)
singleframe_file='finalframes/'+os.path.splitext(basename)[0]+'.coloring_by_structure.gsd'

traj = gsd.open(name=singleframe_file, mode='wb')
traj.append(system)

    


    







