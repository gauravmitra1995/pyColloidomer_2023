import numpy as np
import os
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob
import sys
import argparse

def round_to_int(arr):
    result=[]
    for i in range(arr.size):
        true=arr[i]
        nint=int(arr[i])
        if(np.abs(nint-true)<0.5):
            result.append(np.floor(true))
        else:
            result.append(np.ceil(true))
    return np.array(result)


parser=argparse.ArgumentParser()


#THIS ANALYSIS IS ONLY FOR A TRIMER

#Parameters REQUIRED to be passed as arguments from outside
parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int,default=3)
parser.add_argument("--R",type=float)
parser.add_argument("--Np",type=int)

#Parameters not mandatory to be passed as arguments from outside
parser.add_argument("--metropolis",default=1,type=int)
parser.add_argument("--kAB",type=float,default=200.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--dimension",type=int,default=2)
parser.add_argument("--gammaA",type=float,default=0.1)
parser.add_argument("--gammapatch",type=float,default=0.0001)
parser.add_argument("--r0",type=float,default=2.0)

args = parser.parse_args()
locals().update(vars(args))

if(epsilon!='infinite'):
    epsilon=float(epsilon)

if(Nclusters==2):
    simulationtype='dimer'
elif(Nclusters==3):
    simulationtype='trimer'


print("*******************************************************************************************************************************************************")
print("Np  R  epsilon  r0  gammaA  gammapatch simulationtype")
print(Np," ",R," ",epsilon," ",r0," ",gammaA," ",gammapatch," ",simulationtype)


timesteplist=[]
unbondedbinderslist_allseeds={}
seedlist=[]

timesteps_seeds=[]


#FOR TRIMER, ALSO CALCULATE THE DIFFERENCE IN NUMBER OF UNRECRUITED BINDERS BETWEEN THE TWO TERMINAL DROPLETS AND THE ASYMMETRY PARAMETER

timesteplist=[]
differencelist_allseeds={}
seedlist=[]
meandifferencelist=[]

timesteps_seeds=[]

for filename in sorted(glob.glob(str(os.getcwd())+"/../simulation_setup/"+"/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/selfavoiding/Nclusters"+str(Nclusters)+"/Np"+str(Np)+"/R"+str(R)+"/radiusB"+str(radiusB)+"/kAB"+str(kAB)+"/epsilon"+str(epsilon)+"/dimension"+str(dimension)+"/kspring"+str(kspring)+"/seed*/restlength"+str(r0)+"/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"_seed*.allruns.difference_terminaldroplets.data"),key=lambda x:(int(((os.path.basename(x).split("_")[12]).split(".")[0]).replace('seed','')))):

    seed=int(((os.path.basename(filename).split("_")[12]).split(".")[0]).replace('seed',''))
    data = (np.load(filename,allow_pickle=True))

    unbonded0=np.array(data[1])
    unbonded2=np.array(data[2])

    print(seed,": --->",unbonded0[-1],unbonded2[-1],np.round(data[0][-1],2)," Asymmetry: ",int(data[0][-1]*((unbonded0[-1]+unbonded2[-1])/2.0)))
    
    if(unbonded0[-1]>unbonded0[0] or unbonded2[-1]>unbonded2[0]):
        if unbonded0[-1]!=unbonded2[-1]: 
            #print("Discarding error seed ",seed)
            continue
    
    seedlist.append(seed)
    differencelist_allseeds[seed]=list(np.around(data[0].astype(float),4))
    timesteplist=data[3]
    timesteps_seeds.append(len(timesteplist))

mean_difference_allseeds=[]
stddev_difference_allseeds=[]

print("No of seeds chosen: ",len(seedlist))

if(len(seedlist)>0):
    right_length_array=max(timesteps_seeds)
    temp=[]
    for j in differencelist_allseeds.keys():
        if(len(differencelist_allseeds[j])==right_length_array):
            temp.append(differencelist_allseeds[j])
    mean_difference_allseeds=list(np.mean(temp,axis=0))
    stddev_difference_allseeds=list(np.std(temp,axis=0))
else:
    mean_difference_allseeds=list(np.zeros(data[3].shape[0]))
    stddev_difference_allseeds=list(np.zeros(data[3].shape[0]))

output_file_difference_allseeds=str(os.getcwd())+"/differencepatches_vs_time_data_new/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".averageallseeds.allruns.difference_terminaldroplets.data"

x=np.array([mean_difference_allseeds,stddev_difference_allseeds,timesteplist],dtype=object)
x.dump(output_file_difference_allseeds)



        
        


