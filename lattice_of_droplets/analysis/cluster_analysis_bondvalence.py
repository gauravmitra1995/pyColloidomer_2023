from functions_cluster_analysis import *

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
Np=int((num_particles/num_clusters-1)/2)
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


print('Registering bonded cluster pairs and calculating valence list for each frame...')
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

bond_array_all_frames= np.array(bonds_allframes).astype(int)
timestep_array=np.array(timestep_list,dtype=object)
average_valences=np.around(np.mean(bond_array_all_frames,axis=1),3)


#Save valence vs time data (for all droplets) and averagevalence vs time data 
output_file_bondvalencewithtime = os.path.splitext(trajectory_file)[0]+'.bondvalencewithtime.data'
output_file_averagevalence=os.path.splitext(trajectory_file)[0]+'.averagevalence.data'

r=[]

for i in range(len(timestep_list)):
    p=bond_array_all_frames[i]
    q=timestep_list[i]*dt
    r.append([q,p])

r=np.array(r,dtype=object)
r.dump(output_file_bondvalencewithtime)

y=np.array([timestep_array*dt,average_valences],dtype=object)
y.dump(output_file_averagevalence)

