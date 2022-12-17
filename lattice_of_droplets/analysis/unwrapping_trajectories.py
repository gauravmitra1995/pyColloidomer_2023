from functions_cluster_analysis import *
import hoomd
import hoomd.md

hoomd.context.initialize("")

parser=argparse.ArgumentParser()
parser.add_argument("--trajectory_file",type=str)
parser.add_argument("--minframe",default=0,type=int)
#parser.add_argument("--maxframe",type=int)
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

box_size = np.array(system.configuration.box[:3])


type_list = list(system.particles.types) # list of all particle types in the system
num_particles = system.particles.N # total num of particles in the system
print(num_particles)

A_tags_list=[i for i in range(num_particles) if system.particles.typeid[i]==0] # tags of central particles
num_clusters=len(A_tags_list) # number of clusters
Np=int((num_particles/num_clusters-1)/2)
cluster_ids = np.arange(num_clusters,dtype=int)

typeid_list = [list(system.particles.typeid[i : j]) for i, j in zip([None]+A_tags_list, A_tags_list+[None]) if i!=None] # [[..cluster0 typeids..],[..cluster1 typeids..],..]
Np_list=[(len(i)-1)/2 for i in typeid_list] # list of Np's for each cluster
print(Np_list)

print('Writing unwrapped position (undo PBC) for each non-A particle (only binders) in the system...')

outputfile=os.path.splitext(trajectory_file)[0]+'.unwrap.gsd'
traj=gsd.open(name=outputfile, mode='wb')

for f in range(len(frames_list)):

    system=trajectory[int(frames_list[f])]
    unwrapped_positions=np.zeros((num_particles,3))
    
    for i in range(num_particles):
        center_position=np.array(system.particles.position[i])
        if(system.particles.typeid[i]!=0):
            Aparticle=A_tags_list[find_cluster1(A_tags_list,i)]
            for d in range(3):
                if(center_position[d]-system.particles.position[Aparticle][d]>(box_size[d]/2)):
                    center_position[d]=center_position[d]-box_size[d]
                elif(center_position[d]-system.particles.position[Aparticle][d]<=(-box_size[d]/2)):
                    center_position[d]=center_position[d]+box_size[d]
        unwrapped_positions[i]=np.array(center_position)

    system.particles.position=unwrapped_positions
    traj.append(system)





         



