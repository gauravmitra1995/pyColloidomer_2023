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
systemfin=trajectory[-1]

finaltime=systemfin.configuration.step*dt
box_size = list(system.configuration.box[:3])

type_list = list(system.particles.types) # list of all particle types in the system
num_particles = system.particles.N # total num of particles in the system

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

output_dict1={}
output_dict1['box_size']=box_size
output_dict1['type_list']=type_list
output_dict1['num_particles']=num_particles
output_dict1['A_tags_list']=A_tags_list
output_dict1['num_clusters']=num_clusters
output_dict1['typeid_list']=typeid_list
output_dict1['Np_list']=Np_list
output_dict1['terminal_typeids']=terminal_typeids
output_dict1['num_non_dybonds']=num_non_dybonds
output_dict1['bond_types']=bond_types
output_dict1['dybond_types']=dybond_types
output_dict1['dybond_typeids']=dybond_typeids
output_dict1['dt']=dt
output_dict1['dumpperiod']=trajectory[1].configuration.step-trajectory[0].configuration.step


############################# loop over frames in frame_list to register bonded cluster pairs, num_dybonds #############################
bonded_pairs_allframes=list() # list of lists of bonded cluster pairs in each frame

num_dybonds_allframes=list() # list of number of dybonds in each frame
num_dybonds_by_cluster_allframes=list()


timestep_list=list()

print('Registering bonded cluster pairs and calculating the number of unbonded binders for each droplet.....')
for frame in frames_list:
    system=trajectory[int(frame)]
    timestep_list.append(system.configuration.step) 
    num_dybonds=len(system.bonds.typeid)-num_non_dybonds
    bonded_pairs=list() # list of pairs of bonded clusters for frame frame
    particle_pairs=list()
    if num_dybonds!=0:
        dybond_ids=[i for i in range(len(system.bonds.typeid)) if system.bonds.typeid[i] in dybond_typeids] # bond ids of all dybonds in this frame
        for i in dybond_ids: # iterate over all dybonds (by index)
            pair=system.bonds.group[i] # the the pair of particles (by tags) that form this dybond
            particle1=pair[0] # get the tag of the first particle
            clusterno1=find_cluster1(A_tags_list,particle1) # find which cluster particle 1 belongs to
            particle2=pair[1]
            clusterno2=find_cluster1(A_tags_list,particle2) # find which cluster particle 2 belongs to
            particle_pairs.append([particle1,particle2])
            bonded_pairs.append([clusterno1,clusterno2])
    
    num_dybonds_by_cluster_unsorted=counting_dict(flat_list(bonded_pairs)) # [[1,2],[0,1],[2,1]]-->[1,2,0,1,2,1] --> {1:3, 2:2, 0:1}
    num_dybonds_by_cluster={}
    for i in range(num_clusters):
        if i not in sorted(num_dybonds_by_cluster_unsorted): # if clusterid not in the keys(sorted) of the num_dybonds_by_cluster_unsorted dictionary
            num_dybonds_by_cluster[i]=0
        else:
            num_dybonds_by_cluster[i]=num_dybonds_by_cluster_unsorted[i]
    num_dybonds_by_cluster_allframes.append(num_dybonds_by_cluster)

    bonded_pairs,countrepeats=remove_repeats(bonded_pairs)
    bonded_pairs_allframes.append(bonded_pairs)
    num_dybonds_allframes.append(num_dybonds)


output_dict1['timestep_list']=timestep_list


output_dict2={} # key: frame number, value: a dictionary of number of dybonds
for i in range(len(frames_list)):
    dictionary={}
    dictionary['num_dybonds']=num_dybonds_allframes[i]
    dictionary['num_dybonds_by_cluster']=num_dybonds_by_cluster_allframes[i]
    
    output_dict2[frames_list[i]]=dictionary

output_dict={}
output_dict['general_info']=output_dict1
output_dict['frame_info']=output_dict2


output_filename = os.path.splitext(trajectory_file)[0]+'.npy'
np.save(output_filename, output_dict)

file_path=output_filename

data = (np.load(file_path,allow_pickle=True)).item()

box_size=data['general_info']['box_size']
type_list=data['general_info']['type_list']
num_particles=data['general_info']['num_particles']
A_tags_list=data['general_info']['A_tags_list']
num_clusters=data['general_info']['num_clusters']
typeid_list=data['general_info']['typeid_list']
Np_list=data['general_info']['Np_list']
terminal_typeids=data['general_info']['terminal_typeids']
num_non_dybonds=data['general_info']['num_non_dybonds']
bond_types=data['general_info']['bond_types']
dybond_types=data['general_info']['dybond_types']
dybond_typeids=data['general_info']['dybond_typeids']
timestep_list=data['general_info']['timestep_list']

frames_list=list(data['frame_info'].keys())

num_unbonded_list=list()
num_bonded_list=list()

for i in range(num_clusters): # iterate over all clusters
    num_unbonded_allframes=list()
    num_bonded_allframes=list()
    for frame in frames_list:
        num_unbonded=Np_list[i]-data['frame_info'][frame]['num_dybonds_by_cluster'][i]
        num_bonded=data['frame_info'][frame]['num_dybonds_by_cluster'][i]
        num_unbonded_allframes.append(num_unbonded)
        num_bonded_allframes.append(num_bonded)
    num_unbonded_list.append(num_unbonded_allframes)
    num_bonded_list.append(num_bonded_allframes)

output_num_unbonded = os.path.splitext(trajectory_file)[0]+'.num_unbondedpatches.data'
x=np.array([num_unbonded_list,timestep_list],dtype=object)
x.dump(output_num_unbonded)
