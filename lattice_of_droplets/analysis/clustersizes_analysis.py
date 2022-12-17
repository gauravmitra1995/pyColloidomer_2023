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

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
parser.add_argument("--minframe",default=0,type=int)
#parser.add_argument("--maxframe",default=100,type=int)
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

bonds_allframes=[]

A_tags_list=[i for i in range(num_particles) if system.particles.typeid[i]==0] # tags of central particles
num_clusters=len(A_tags_list) # number of clusters
#Np=int((num_particles/num_clusters-1)/2)
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

timestep_list=[]
time_list=[]

print('Registering bonded cluster pairs and performing cluster analysis of various structures for each frame...')

frames_list=frames_list[:]

result_df=pd.DataFrame()

structures=['monomers','dimers','linear chains','loops','other']
results_list=[]
linearchainlengths_allframes=[]

for frame in frames_list:
    fraction_structures={}
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
    for no in range(len(cluster_ids)):
        id=cluster_ids[no]
        for x in range(len(cluster_bondtable_new)):
            clusterno1=cluster_bondtable_new[x][0]
            clusterno2=cluster_bondtable_new[x][1]
            if(id==clusterno1 or id==clusterno2):
                 number_bonds[id]=number_bonds[id]+1
    bonds_allframes.append(number_bonds)
    timestep_list.append(system.configuration.step)
    time_list.append(system.configuration.step*dt)

    number_of_droplets_as_monomers=len(monomers)
    number_of_droplets_in_linear_chains_not_dimers=0
    number_of_droplets_in_dimers=0
    number_of_droplets_in_loops=0
    number_of_droplets_in_gelsorbranchedchains=0

    linear_chain_sizes=[]

    for cluster in clustering:
        avg_valence=0
        for i in cluster:
            avg_valence=avg_valence+int(number_bonds[int(i)])
        avg_valence=avg_valence/len(cluster)
        avg_valence=np.round(avg_valence,3)
        avg_valence_linear_chain=np.round((1.0*2 + (len(cluster)-2)*2.0)/len(cluster),3)
       
        valences=[]
        for i in cluster:
            valence_i=int(number_bonds[int(i)])
            valences.append(valence_i)

        no_of_1s=countX(valences,1)
        no_of_2s=countX(valences,2)

        if(len(cluster)==2):
            if(no_of_1s==len(cluster)):
                number_of_droplets_in_dimers+=len(cluster)
        elif(len(cluster)>2):
            if(no_of_1s==0 and no_of_2s==len(cluster)):
                number_of_droplets_in_loops+=len(cluster)
            elif(no_of_1s==2 and no_of_2s==len(cluster)-2):
                number_of_droplets_in_linear_chains_not_dimers+=len(cluster)
                linear_chain_sizes.append(len(cluster))
            else:
                number_of_droplets_in_gelsorbranchedchains+=len(cluster)
                
    if(len(linear_chain_sizes)!=0):
        linearchainlengths_allframes.append(np.array(linear_chain_sizes,dtype=object))


    fraction_structures['monomers']=np.round(number_of_droplets_as_monomers/len(A_tags_list),3)
    fraction_structures['dimers']=np.round(number_of_droplets_in_dimers/len(A_tags_list),3)
    fraction_structures['linear chains']=np.round(number_of_droplets_in_linear_chains_not_dimers/len(A_tags_list),3)
    fraction_structures['loops']=np.round(number_of_droplets_in_loops/len(A_tags_list),3)
    fraction_structures['other']=np.round(number_of_droplets_in_gelsorbranchedchains/len(A_tags_list),3)

    results_list.append(fraction_structures)
    
result_df=pd.DataFrame(results_list)
print(result_df)

result_df.to_csv(os.path.splitext(trajectory_file)[0]+".structuralanalysis.csv",sep=" ",index=False)

linearchainlengths_allframes=np.array(linearchainlengths_allframes,dtype=object)
output_file_linearchainlengths = os.path.splitext(trajectory_file)[0]+'.linearchainlengths.data'
linearchainlengths_allframes.dump(output_file_linearchainlengths)

