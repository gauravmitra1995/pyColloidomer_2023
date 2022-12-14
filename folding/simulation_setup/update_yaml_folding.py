import numpy as np
import sys
import os
import os.path
import argparse
import yaml
import copy


with open ('input_general.yaml') as f:
        data_general = yaml.load(f,Loader=yaml.FullLoader)

with open ('input_clusters.yaml') as f:
        data_clusters = yaml.load(f,Loader=yaml.FullLoader)

with open ('input_particles.yaml') as f:
        data_particles = yaml.load(f,Loader=yaml.FullLoader)

parser = argparse.ArgumentParser()
parser.add_argument('--Np', action= "store",nargs ="*", default=500, dest = "Np_per_sup", type = int)
parser.add_argument('--R', action= "store",nargs ="*", default=50.0, dest = "R", type = float)
parser.add_argument('--simulationtype', action= "store", nargs="*", default='polymer', dest="simulationtype", type=str)
parser.add_argument('--radiusB',action="store", nargs="*", default=1.0, dest="radiusB",type=float)
parser.add_argument('--Nclusters',action="store", nargs="*", default=2, dest="Nclusters",type=int)
parser.add_argument('--seed',action="store", nargs="*", default=10, dest="seed",type=int)
parser.add_argument('--koninit',action="store", nargs="*", default=100.0, dest="koninit",type=float)
parser.add_argument('--koffinit',action="store", nargs="*", default=0, dest="koffinit",type=float)
parser.add_argument('--dimension',action="store", nargs="*", default=2, dest="dimension", type=int)
parser.add_argument('--kAB',action="store", nargs="*", default=200.0, dest="kAB", type=float)
parser.add_argument('--gammaA',action="store", nargs="*", default=0.1, dest="gammaA",type=float)
parser.add_argument('--gammapatch',action="store", nargs="*", default=0.0001, dest="gammapatch",type=float)
parser.add_argument('--metropolis',action="store", nargs="*",default=1,dest="metropolis",type=int)
parser.add_argument('--kspring',action="store", nargs="*", default=10.0, dest="kspring", type=float)
parser.add_argument('--zerodynbondlength', action= "store", nargs="*", default='False', dest="zerodynbondlength", type=str)
parser.add_argument('--kT',action= "store", nargs="*",default=1.0,dest="kT",type=float)


print("")
print("Adjusted Values")
print("--------------------")
results = parser.parse_args()
print("Simulationtype:",results.simulationtype)
print("Np_per_sup:",results.Np_per_sup)
print("R:",results.R)
print("radiusB:",results.radiusB)
print("Nclusters:",results.Nclusters)
print("seed:",results.seed)
print("koninit:",results.koninit)
print("koffinit:",results.koffinit)
print("dimension:",results.dimension)
print("AB spring constant:",results.kAB)
print("gammaA:",results.gammaA)
print("gammapatch:",results.gammapatch)
print("Metropolis flag:",results.metropolis)
print("Dynamic bond spring constant:",results.kspring)
print("Zero bond length flag:",results.zerodynbondlength)
print("kT:",results.kT)

simulationtype=results.simulationtype[0]

result_list = []

keys_general=list(data_general.keys())
keys_clusters=list(data_clusters.keys())
keys_particles=list(data_particles.keys())


for Np_per_sup in results.Np_per_sup:
    for R in results.R:
        for radiusB in results.radiusB: 
            for Nclusters in results.Nclusters:
                for seed in results.seed:
                    for koninit in results.koninit:
                        for koffinit in results.koffinit:
                            for dimension in results.dimension:
                                for kAB in results.kAB:
                                    for gammaA in results.gammaA:
                                        for gammapatch in results.gammapatch:
                                            for metropolis in results.metropolis:
                                                for kspring in results.kspring:
                                                    for zerodynbondlength in results.zerodynbondlength:
                                                        for kT in results.kT:

                                                            new_data_clusters = copy.copy(data_clusters)
                                                            new_data_particles = copy.copy(data_particles)
                                                            new_data_general = copy.copy(data_general)        

                                                     
                                                            if(koffinit==0.0):
                                                                epsilon='infinite'
                                                            else:
                                                                if(koninit!=0.0):
                                                                    epsilon=format(np.log(koninit/koffinit),'.1f')
                                                                else:
                                                                    epsilon=0.0

                                                            if(zerodynbondlength=='True'):
                                                                restlength=0.0
                                                                CCsoftVepsilon=0
                                                            else:
                                                                restlength=2.0*radiusB
                                                                CCsoftVepsilon=200.0

                                                            
                                                            for x in keys_clusters: 
                                                                new_data_clusters[x]['Nc']=int(Nclusters)
                                                                new_data_clusters[x]['Np_per_sup']= Np_per_sup


                                                            for y in keys_particles:
                                                                if(y=='particle'):
                                                                    new_data_particles[y]['A']['radius']=R
                                                                    new_data_particles[y]['B']['radius']=radiusB
                                                                if(y=='soft_V'):
                                                                    new_data_particles[y]['CC']['epsilon']=CCsoftVepsilon
                                                                    new_data_particles[y]['DD']['epsilon']=CCsoftVepsilon

                                                            for y in keys_particles:
                                                                if(y=='harmonic'):
                                                                    new_data_particles[y]['AB']['kspring']=kAB

                                                            for y in keys_particles:
                                                                if(y=='dybond'):
                                                                    for i in data_particles[y].keys():
                                                                        bond_dict=data_particles[y][i]
                                                                        key=i
                                                                        new_data_particles[y][key]['kon_init']=koninit
                                                                        if(key=='DD'):
                                                                            new_data_particles[y][key]['koff_init']=koffinit
                                                                            new_data_particles[y][key]['self_avoid_chain']=0
                                                                        elif(key=='CC'):
                                                                            new_data_particles[y][key]['koff_init']=0.0
                                                                            new_data_particles[y][key]['self_avoid_chain']=1

                                                                        new_data_particles[y][key]['metropolis']=metropolis
                                                                        new_data_particles[y][key]['kspring']=kspring
                                                                                
                                                            new_data_general['simulationtype']=simulationtype
                                                            new_data_general['seed']=seed
                                                            new_data_general['dimension']=dimension
                                                            new_data_general['gammaA']=gammaA
                                                            new_data_general['gammapatch']=gammapatch
                                                            new_data_general['lattice_dim']=[1,int(Nclusters)]
                                                            new_data_general['zerodynbondlength']=zerodynbondlength
                                                            new_data_general['kT']=kT


                                                            new_data_general['fileprefix']=str(simulationtype)+"_restl"+str(restlength)+"_Nc"+str(Nclusters)+"_Np"+str(Np_per_sup)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"_seed"+str(seed)

                                                            path=os.path.join(os.getcwd(),str(simulationtype),"gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch),"Nclusters"+str(Nclusters),"Np"+str(Np_per_sup),"R"+str(R),"radiusB"+str(radiusB),"kAB"+str(kAB),"epsilon"+str(epsilon),"dimension"+str(dimension),"kspring"+str(kspring),"seed"+str(seed),"restlength"+str(restlength))

                                                            file_general=os.path.join(path,"input_general.yaml")
                                                            file_clusters=os.path.join(path,"input_clusters.yaml")
                                                            file_particles=os.path.join(path,"input_particles.yaml")
                                        
                                                            os.makedirs(os.path.dirname(file_general), exist_ok=True)
                                                            os.makedirs(os.path.dirname(file_clusters), exist_ok=True)
                                                            os.makedirs(os.path.dirname(file_particles), exist_ok=True) 
                                        
                                                            with open (file_general, 'w') as outfile_general:
                                                                yaml.dump(new_data_general, outfile_general, default_flow_style= False)
                                                            print("Updated Yaml General input File generated.")

                                                            with open (file_clusters, 'w') as outfile_clusters:
                                                                yaml.dump(new_data_clusters, outfile_clusters, default_flow_style= False)
                                                            print("Updated Yaml Cluster input File generated.")

                                                            with open (file_particles, 'w') as outfile_particles:
                                                                yaml.dump(new_data_particles, outfile_particles, default_flow_style= False)
                                                            print("Updated Yaml Particle input File generated.")
                                        
                                                            print("*******************************************************************************************")
                                        
                                                            print("Full path to run directory:")
                                                            print(os.path.dirname(file_general))

