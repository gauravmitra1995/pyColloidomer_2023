from functions_cluster_analysis import *

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
parser.add_argument("--minframe",default=0,type=int)
#parser.add_argument("--maxframe",default=100,type=int)
parser.add_argument("--dt",default=0.0005,type=float)
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

print(A_tags_list)

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

mean_noofbonds_adjacentpairs_allframes_CC=[]
mean_noofbonds_adjacentpairs_allframes_DD=[]

num_DDbonds_allframes=[]
num_CCbonds_allframes=[]

print('Registering number of binders in the adhesion patches between adjacent droplets only (for C and D types) ...')

for frame in frames_list[:]:
    system=trajectory[int(frame)]
    bondtable=system.bonds.group.astype(int)
    cluster_bondtable=[]
    dybondtable=[]
    cluster_bondtable_dybonds=[]

    for i in range(bondtable.shape[0]):
        bond=bondtable[i]
        particle1=bond[0]
        particle2=bond[1]
        cluster1=int(find_cluster1(A_tags_list,particle1))
        cluster2=int(find_cluster1(A_tags_list,particle2))
        
        clusterpair=[cluster1,cluster2]
        cluster_bondtable.append(clusterpair)
        
        if(clusterpair[0]!=clusterpair[1]):    #filter out pairs where the two droplet id's in a pair are not same
            dybondtable.append(bondtable[i])
            cluster_bondtable_dybonds.append(clusterpair)

    dybondtable=np.array(dybondtable,dtype=object)
    cluster_bondtable_dybonds_unique,countrepeats=remove_repeats(cluster_bondtable_dybonds)    #remove repeats to get unique pairs of bonded droplets 

    cluster_bondtable_dybonds_unique_onlyadjacent=[]

    for c in cluster_bondtable_dybonds_unique:
        pair=c
        if(np.abs(pair[1]-pair[0])==1): #check if these droplets under consideration are adjacent or not 
            cluster_bondtable_dybonds_unique_onlyadjacent.append(c)

    type_table=[]
    dybondtable_onlyDD=[]
    dybondtable_onlyCC=[]
    clusterbondtable_onlyDD=[]
    clusterbondtable_onlyCC=[]
    for d in dybondtable:
        p1=d[0]
        p2=d[1]

        #get bond table of droplet pairs for C-C and D-D separately 
        if(system.particles.typeid[p1]==2 and system.particles.typeid[p2]==2):
            type_table.append(['C','C'])
            dybondtable_onlyCC.append(list(d))
            clus1=int(find_cluster1(A_tags_list,p1))
            clus2=int(find_cluster1(A_tags_list,p2))
            clusterbondtable_onlyCC.append([clus1,clus2])
        elif(system.particles.typeid[p1]==3 and system.particles.typeid[p2]==3):
            type_table.append(['D','D'])
            dybondtable_onlyDD.append(list(d))
            clus1=int(find_cluster1(A_tags_list,p1))
            clus2=int(find_cluster1(A_tags_list,p2))
            clusterbondtable_onlyDD.append([clus1,clus2])

    num_DDbonds_allframes.append(len(dybondtable_onlyDD))  #no of D-D bonds 
    num_CCbonds_allframes.append(len(dybondtable_onlyCC))  #no of C-C bonds

    #Count no of adjacent C-C and D-D bonds 
    count_DD_adjacent=0
    count_CC_adjacent=0
    for c in cluster_bondtable_dybonds_unique_onlyadjacent:
        for cbt in clusterbondtable_onlyCC:
            if(c==[cbt[0],cbt[1]] or c==[cbt[1],cbt[0]]):
                count_CC_adjacent+=1
        for dbt in clusterbondtable_onlyDD:
            if(c==[dbt[0],dbt[1]] or c==[dbt[1],dbt[0]]):
                count_DD_adjacent+=1

    #to find mean no of adjacent C-C  / D-D bonds : divide the total by the number of adhesion patches between adjacent pairs of droplets (in this case 6, for 7 droplets)

    mean_noofbonds_adjacentpairs_allframes_CC.append(count_CC_adjacent/len(cluster_bondtable_dybonds_unique_onlyadjacent))
    mean_noofbonds_adjacentpairs_allframes_DD.append(count_DD_adjacent/len(cluster_bondtable_dybonds_unique_onlyadjacent))

    timestep_list.append(system.configuration.step)


print("*************************************************************************************************************************************************************************************************************")
print("Mean no of bonds in adhesion patch (for a given bond type) for a droplet pair:")
print(mean_noofbonds_adjacentpairs_allframes_CC)
print(mean_noofbonds_adjacentpairs_allframes_DD)

time=np.array(timestep_list)*dt
l=time.shape[0]

fig,ax=plt.subplots(figsize=(20,15),dpi=100)
ax.plot(time[int(l/2):][::3],mean_noofbonds_adjacentpairs_allframes_CC[int(l/2):][::3],marker='.',markersize=15.0,linewidth=5.0,color='red',label='C-C')
ax.plot(time[int(l/2):][::3],mean_noofbonds_adjacentpairs_allframes_DD[int(l/2):][::3],marker='.',markersize=15.0,linewidth=5.0,color='darkblue',label='D-D')
ax.legend(loc='upper right',prop={'size':40},ncol=2)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(5)
ax.tick_params(labelsize=50,axis='y',pad=20)
ax.tick_params(labelsize=40,axis='x',pad=20)
ax.set_xticks(np.arange(80000,150000,20000))
ax.set_yticks(np.arange(0,110,10))
ax.set_ylabel(r'${\langle N_{\mathrm{bonds}} \rangle}^{\mathrm{adj}}$',labelpad=20,fontsize=70)
ax.set_xlabel('Simulation time (in HOOMD units)',labelpad=20,fontsize=60)
fig.tight_layout()
plt.savefig('final_figures/num_CC_and_DD_bonds.png',bbox_inches='tight')
plt.close()

