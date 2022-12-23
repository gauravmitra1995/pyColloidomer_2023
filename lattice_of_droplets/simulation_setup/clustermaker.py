import numpy as np
class patchyparticlecluster:

    def __init__(self,Np_per_sup,typelist,centerpos,clusterid,r_dict,mass_dict,input_r0_dict,choice_initialdistribution):

            #initializes the attributes of each object (droplet with patches on surface)

            self.Np_per_sup=Np_per_sup # number of patches (B-X pairs) per superpatch
            self.Nsup=1 #no of superpatches in each cluster (right now, it is 1 for all)
            self.typelist=typelist
            self.R=r_dict[typelist[0]] # Assuming the first element of typelist is the central A particle
            self.clusterid=clusterid

            self.sigmaB=r_dict[typelist[1]]*2
            self.massB=mass_dict[typelist[1]]
            self.r0AB=input_r0_dict[''.join([self.typelist[0][0],self.typelist[1][0]])] # ['A','B']--> 'AB'

            self.sigmaX=r_dict[typelist[2]]*2
            self.massX=mass_dict[typelist[2]]
            self.r0AX=input_r0_dict[''.join([self.typelist[0][0],self.typelist[2][0]])] # ['A','X']--> 'AX'

            self.choice_initialdistribution=choice_initialdistribution
            self.phimax=np.pi # opening angle of each superpatch (0 to pi)
            self.Np=self.Np_per_sup*self.Nsup # total number of patches (B-X pairs) in the cluster
            self.vmin=(np.cos(self.phimax)+1)/2 # corresponds to phimax, for assigning patches randomly
            self.N=self.Np*2+1 # total number of particles in the cluster

            self.particle_positions=[centerpos]+[None]*self.Np*2
            self.particle_types=[typelist[0]]+[None]*self.Np*2 #1st element in typelist is central particle (typelist[0])
            self.particle_diameters=[self.R*2]+[None]*self.Np*2
            self.particle_masses=[mass_dict[typelist[0]]]+[None]*self.Np*2
            self.particle_clusterids=[self.clusterid]+[None]*self.Np*2
            self.particle_typeids=[0]+[None]*self.Np*2

            self.centerpos=centerpos
            self.bonds={}
            self.angles={}



    def assign_patches(self,chainlink):

        #outputs coordinates, diameters, and types of the patches and store in object attributes

        import random
        from read_parameters import spherical2cart,sphere_fibonacci_grid_points

        if(self.choice_initialdistribution==2):
            pointsB,pointsX= sphere_fibonacci_grid_points ( self.Np,(self.r0AB),(self.r0AX) )

        for i in range(1,self.Np+1):
            u=random.uniform(0,1)
            theta=2*np.pi*u
            v=random.uniform(self.vmin,1)
            phi=np.arccos(2*v-1)

            if self.Nsup==1:
                if(self.choice_initialdistribution==1):
                    x1,y1,z1=spherical2cart((self.r0AB),theta,phi)
                    x2,y2,z2=spherical2cart((self.r0AX),theta,phi)
                elif(self.choice_initialdistribution==2):
                    x1,y1,z1=pointsB[i-1][0],pointsB[i-1][1],pointsB[i-1][2]
                    x2,y2,z2=pointsX[i-1][0],pointsX[i-1][1],pointsX[i-1][2]


            self.particle_positions[i]=list(np.add([x1,y1,z1],self.centerpos))
            self.particle_types[i]=self.typelist[1]
            self.particle_diameters[i]=self.sigmaB
            self.particle_masses[i]=self.massB
            self.particle_clusterids[i]=self.clusterid

            self.particle_positions[i+self.Np]=list(np.add([x2,y2,z2],self.centerpos))
            self.particle_types[i+self.Np]=self.typelist[2]
            self.particle_diameters[i+self.Np]=self.sigmaX
            self.particle_masses[i+self.Np]=self.massX
            self.particle_clusterids[i+self.Np]=self.clusterid

        idx=len(self.typelist)-1
        d=len(self.typelist)-2
        ctr=0
        while(ctr<(len(self.typelist)-3)):
            for i in range(int(idx*self.Np/d)+1,int((idx+1)*self.Np/d)+1):
                self.particle_types[i]=self.typelist[3+ctr]
            idx+=1
            ctr+=1

        if chainlink =='True' and self.Np > 1:

            from read_parameters import rotation

            for i in range(1,3): # replace the 1st and 2nd patches with one on left and one on right

                if i==1:  # place a patch pair (B and X) on the left
                    Mleft=rotation((np.pi/2),'x')
                    x1,y1,z1=Mleft.dot(spherical2cart((self.r0AB),0,0))
                    x2,y2,z2=Mleft.dot(spherical2cart((self.r0AX),0,0))

                elif i==2: # place a patch pair (B and X) on the right
                    Mright=rotation((-np.pi/2),'x')
                    x1,y1,z1=Mright.dot(spherical2cart((self.r0AB),0,0))
                    x2,y2,z2=Mright.dot(spherical2cart((self.r0AX),0,0))

                self.particle_positions[i]=list(np.add([x1,y1,z1],self.centerpos))
                self.particle_types[i]=self.typelist[1]
                self.particle_diameters[i]=self.sigmaB
                self.particle_masses[i]=self.massB
                self.particle_clusterids[i]=self.clusterid


                self.particle_positions[i+self.Np]=list(np.add([x2,y2,z2],self.centerpos))
                self.particle_types[i+self.Np]=self.typelist[2]
                self.particle_diameters[i+self.Np]=self.sigmaX
                self.particle_masses[i+self.Np]=self.massX
                self.particle_clusterids[i+self.Np]=self.clusterid

    def set_bonds_angles(self,previousN):
        """
        creates a dictionary of bonds and angles, where key: bond pair (e.g. 'AB') and value: list of bonding pairs
        'previousN': e.g. if the first cluster has N=61 (2*Np+1, where Np=30) and we are calling set_bonds_angles for the second cluster, then the IDs/tags
        of the particles in this second cluster are (61,62,...)
        """
        self.bonds[self.typelist[0][0]+self.typelist[1][0]] = [[previousN+0,previousN+i] for i in range(1,self.Np+1)] #defining A-B Bonds
        l=len(self.typelist)
        d=l-2
        idx=0
        while(idx<=(l-3)):
            self.bonds[self.typelist[1][0]+self.typelist[idx+2][0]] = [[previousN+i,previousN+i+self.Np] for i in range(int(idx*self.Np/d)+1,int((idx+1)*self.Np/d)+1)] #defining B-X Bonds
            self.angles[self.typelist[0][0]+self.typelist[1][0]+self.typelist[idx+2][0]] = [[previousN+0,previousN+i,previousN+i+self.Np] for i in range(int(idx*self.Np/d)+1,int((idx+1)*self.Np/d)+1)]
            #defining A-B-X angles
            idx+=1

def main():
    pass

if __name__ == '__main__':
    main()

