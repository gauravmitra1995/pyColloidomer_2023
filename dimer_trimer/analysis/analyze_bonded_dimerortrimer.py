from gsd import hoomd as gsd
import numpy as np
import sys

def get_dynbond_trj(trj):
    dynbond_list = []
    for snapshot in trj:
        dynbond_list.append(get_dynbonded(snapshot))
    return dynbond_list 

def get_dynbonded(snapshot):
    bonds = snapshot.bonds.group
    dynbond_type = len(snapshot.bonds.types)-1
    dynbonds = bonds[snapshot.bonds.typeid==dynbond_type,:]
    dynbonds = np.sort(dynbonds,axis=1)
    return dynbonds

def align_xyz_vec(xyz, droplet_vec, target_vec=np.array((0,1,0))):
    #https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    droplet_vec = droplet_vec/np.linalg.norm(droplet_vec)

    v = np.cross(droplet_vec, target_vec)
    c = np.dot(droplet_vec, target_vec)

    h = (1 - c)/(1 - c**2)

    vx, vy, vz = v
    rot =[[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
          [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
          [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]]
    xyz = np.dot(rot,xyz.T).T

    return xyz

def align_trj_droplet(trj,droplet_type='A'):
    for i in range(len(trj)):
        snapshot = trj[i]
        droplet_typeid = snapshot.particles.types.index(droplet_type)
        droplet_index = np.where(snapshot.particles.typeid==droplet_typeid)[0]
        droplet_positions = snapshot.particles.position[droplet_index]
        droplet_vec = droplet_positions[1] - droplet_positions[0]
        new_xyz = align_xyz_vec(snapshot.particles.position, droplet_vec)
        new_droplet_positions = new_xyz[droplet_index]
        new_com = new_droplet_positions.mean(axis=0)
        new_xyz -= new_com
        snapshot.particles.position = new_xyz
    return trj


def dynbond_gsd(gsd_file,droplet_type='A'):
    try: 
        gsd_fh = gsd.open(gsd_file,'rb')
    except RuntimeError:
        return None,None,None
    trj = [frame for frame in gsd_fh]
    dynbond_trj = get_dynbond_trj(trj)
    trj = align_trj_droplet(trj,droplet_type=droplet_type)
    droplet_typeid = trj[0].particles.types.index(droplet_type)
    droplet_ids = np.where(trj[0].particles.typeid==droplet_typeid)[0]
    #this is for a dimer
    Np = int((droplet_ids[1]-1)/2)

    return dynbond_trj, Np, trj

