import numpy as np
import sys
import os
import os.path
import argparse
import yaml
import copy

parser = argparse.ArgumentParser()
parser.add_argument('--Np', action= "store",nargs ="*", default=500, dest = "Np_per_sup", type = int)
parser.add_argument('--R', action= "store",nargs ="*", default=100.0, dest = "R", type = float)
parser.add_argument('--radiusB',action="store", nargs="*", default=1.0, dest="radiusB",type=float)
parser.add_argument('--simulationtype', action= "store", nargs="*", default='lattice', dest="simulationtype", type=str)
parser.add_argument('--Nclusters',action="store", nargs="*", default=64, dest="Nclusters",type=int)
parser.add_argument('--areafraction',action="store", nargs="*", default=0.4, dest="areafraction",type=float)
parser.add_argument('--seed',action="store", nargs="*", default=10, dest="seed",type=int)
parser.add_argument('--koninit',action="store", nargs="*", default=100.0, dest="koninit",type=float)
parser.add_argument('--koffinit',action="store", nargs="*", default=0, dest="koffinit",type=float)
parser.add_argument('--gammaA',action="store", nargs="*", default=0.1, dest="gammaA",type=float)
parser.add_argument('--gammapatch',action="store", nargs="*", default=0.0001, dest="gammapatch",type=float)
parser.add_argument('--metropolis',action="store", nargs="*",default=1,dest="metropolis",type=int)
parser.add_argument('--kspring',action="store",nargs="*",default=10.0,dest="kspring",type=float)
parser.add_argument('--zerodynbondlength', action= "store", nargs="*", default='False', dest="zerodynbondlength", type=str)
print("")
print("Adjusted Values")
print("--------------------")
results = parser.parse_args()
print("Np_per_sup:",results.Np_per_sup)
print("R:",results.R)
print("radiusB:",results.radiusB)
print("Nclusters:",results.Nclusters)
print("areafraction:",results.areafraction)
print("seed:",results.seed)
print("koninit:",results.koninit)
print("koffinit:",results.koffinit)
print("gammaA:",results.gammaA)
print("gammapatch:",results.gammapatch)
print("Metropolis flag:",results.metropolis)
print("Dynbond kspring:",results.kspring)
print("Zero bond length flag:",results.zerodynbondlength)

result_list = []
simulationtype=results.simulationtype[0]


if(simulationtype=='lattice'):
    with open ('input_general_lattice.yaml') as f:
            data_general = yaml.load(f,Loader=yaml.FullLoader)
    with open ('input_clusters_lattice.yaml') as f:
            data_clusters = yaml.load(f,Loader=yaml.FullLoader)
else:
    with open ('input_general_CN.yaml') as f:
            data_general = yaml.load(f,Loader=yaml.FullLoader)
    with open ('input_clusters_CN.yaml') as f:
            data_clusters = yaml.load(f,Loader=yaml.FullLoader)

with open ('input_particles.yaml') as f:
        data_particles = yaml.load(f,Loader=yaml.FullLoader)

keys_general=list(data_general.keys())
keys_clusters=list(data_clusters.keys())
keys_particles=list(data_particles.keys())

for Np_per_sup in results.Np_per_sup:
    for R in results.R:
        for radiusB in results.radiusB: 
            for Nclusters in results.Nclusters: 
                for areafraction in results.areafraction:
                    for seed in results.seed:
                        for koninit in results.koninit:
                            for koffinit in results.koffinit:
                                for gammaA in results.gammaA:
                                    for gammapatch in results.gammapatch:
                                        for metropolis in results.metropolis:
                                            for kspring in results.kspring:
                                                for zerodynbondlength in results.zerodynbondlength:

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
                                                        CDsoftVepsilon=0
                                                    else:
                                                        restlength=2.0*radiusB
                                                        CDsoftVepsilon=200.0

                                                    
                                                    idx=0
                                                    for x in keys_clusters: 
                                                        if(simulationtype=='coordination_no'):
                                                            if(new_data_clusters[x]['Np_per_sup'] is None):
                                                                new_data_clusters[x]['Np_per_sup']= Np_per_sup
                                                        else:
                                                            new_data_clusters[x]['Np_per_sup']= Np_per_sup
                        
                                                        if(type(new_data_clusters[x]['Nc']) is not list):
                                                            if(idx==0):
                                                                new_data_clusters[x]['Nc']=int(Nclusters/2)
                                                            elif(idx==1):
                                                                if(Nclusters%2==0):
                                                                    new_data_clusters[x]['Nc']=int(Nclusters/2)
                                                                else:
                                                                    new_data_clusters[x]['Nc']=int(Nclusters/2) + 1
                                                        if(idx==1):
                                                            break
                                                        idx+=1

                                                    for y in keys_particles:
                                                        if(y=='particle'):
                                                            new_data_particles[y]['A']['radius']=R
                                                            new_data_particles[y]['B']['radius']=radiusB
                                                        if(y=='soft_V'):
                                                            new_data_particles[y]['CD']['epsilon']=CDsoftVepsilon
                                                   

                                                    for y in keys_particles:
                                                        if(y=='dybond'):
                                                            for i in data_particles[y].keys():
                                                                bond_dict=data_particles[y][i]
                                                                key=i  
                                                                new_data_particles[y][key]['kon_init']=koninit
                                                                new_data_particles[y][key]['koff_init']=koffinit
                                                                new_data_particles[y][key]['metropolis']=metropolis
                                                                new_data_particles[y][key]['kspring']=kspring
                                                                    
                                                    new_data_general['simulationtype']=simulationtype
                                                    new_data_general['areafraction']=areafraction
                                                    new_data_general['seed']=seed
                                                    new_data_general['gammaA']=gammaA
                                                    new_data_general['gammapatch']=gammapatch
                                                    new_data_general['zerodynbondlength']=zerodynbondlength

                                                    if(simulationtype=='lattice'):
                                                        new_data_general['lattice_dim']=[int(np.sqrt(Nclusters)),int(np.sqrt(Nclusters))]

                                                    Nclusters=new_data_general['lattice_dim'][0]*new_data_general['lattice_dim'][1]
                                                    Np_per_sup=new_data_clusters['ABD']['Np_per_sup']

                                                    new_data_general['fileprefix']=str(simulationtype)+"_restl"+str(restlength)+"_Nc"+str(Nclusters)+"_Np"+str(Np_per_sup)+"_R"+str(R)+"_rB"+str(radiusB)+"_areafrac"+str(areafraction)+"_eps"+str(epsilon)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"_seed"+str(seed)
                                                    path=os.path.join(os.getcwd(),str(simulationtype),"gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch),"Nclusters"+str(Nclusters),"Np"+str(Np_per_sup),"R"+str(R),"radiusB"+str(radiusB),"areafraction"+str(areafraction),"epsilon"+str(epsilon),"kspring"+str(kspring),"seed"+str(seed),"restlength"+str(restlength))
                                                      
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

