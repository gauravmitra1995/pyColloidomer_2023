import numpy as np
import os
import pickle
from gsd import hoomd as gsd
import matplotlib.pyplot as plt
import glob
import sys
import seaborn as sns
sns.set_context("poster", font_scale=1.0)
import argparse
import pandas as pd

def round_to_int(val):

    true=val
    nint=int(val)
    if(np.abs(nint-true)<0.5):
        rounded=np.floor(true)
    else:
        rounded=np.ceil(true)
    return np.array(rounded)

parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside
parser.add_argument("--epsilon",default=20.7)
parser.add_argument("--Nclusters",type=int)
#parser.add_argument("--R",type=float)
#parser.add_argument("--Np",type=int)

#Parameters not mandatory to be passed as arguments from outside
parser.add_argument("--metropolis",default=1,type=int)
parser.add_argument("--kAB",type=float,default=200.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--dimension",type=int,default=2)
parser.add_argument("--gammaA",type=float,default=0.1)
parser.add_argument("--gammapatch",type=float,default=0.0001)
parser.add_argument("--r0",type=float,default=2.0)
parser.add_argument("--dt",default=0.001,type=float)

args = parser.parse_args()
locals().update(vars(args))

if(epsilon!='infinite'):
    epsilon=float(epsilon)

def get_phasediagrams(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid):

    Nplist=[]
    Rlist=[]
    saturationfractionunbonded_dict={}
    saturationfractionbonded_dict={}
    results=[]
    
    for filename in sorted(glob.glob(str(os.getcwd())+"/saturation_data/saturationfraction_restl"+str(r0)+"_Nc"+str(Nclusters)+"_clusterid"+str(clusterid)+"_Np*_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(epsilon)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".varyingRandNp.data"),key=lambda x:int(os.path.basename(x).split("_")[4].replace('Np',''))): 
        data=np.load(filename,allow_pickle=True)
        Np=int(os.path.basename(filename).split("_")[4].replace('Np',''))
        if(Np not in Nplist):
            Nplist.append(Np)
        Rlist=data[0]
        saturationfractionlist=[]
        saturationfractionlist=data[1]
        
        for i in range(len(Rlist)):
            R=Rlist[i]
            label=(R,Np)
            saturationfrac=saturationfractionlist[i]
            saturationfractionunbonded_dict[label]=saturationfrac
            saturationfractionbonded_dict[label]=1.0-saturationfrac

    for i in saturationfractionunbonded_dict.keys():
        label=i
        Np=label[1]
        R=label[0]
        saturationfrac=saturationfractionunbonded_dict[label]
        results.append(np.array([R,Np,saturationfrac]))

    results=np.array(results)
    """    
    ctr=0
    fig,ax=plt.subplots(figsize=(20,15),dpi=100)
    colormap=plt.cm.jet
    colors = [colormap(i) for i in np.linspace(0,1,len(Rlist))]
    r=0
    for R in Rlist:
        saturationfracbonded_allNp=[]
        saturationfracunbonded_allNp=[]
        numbonded_allR=[]
        for Np in Nplist:
            label=(R,Np)
            saturationfracunbonded=saturationfractionunbonded_dict[label]
            saturationfracbonded=1.0-saturationfracunbonded
            numbonded=round_to_int(saturationfracbonded*label[0])
            saturationfracbonded_allNp.append(saturationfracbonded)
            saturationfracunbonded_allNp.append(saturationfracunbonded)
            numbonded_allR.append(numbonded)
            ctr+=1
        plt.plot(Nplist,saturationfracunbonded_allNp,marker='o',markersize=20.0,color=colors[r],linewidth=10.0,label=r'$R=$'+str(R))
        print(saturationfracunbonded_allNp)
        print(Nplist)
        r+=1
    
    plt.legend(loc='best',ncol=2,prop={'size': 40})
    plt.xlabel(r'$N_{b}$',fontsize=60)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    if(Nclusters==2 and clusterid==0):
        plt.title('Dimer',fontsize=60)
    elif(Nclusters==2 and clusterid==1):
        plt.title('Dimer',fontsize=60)
    elif(Nclusters==3 and clusterid==0):
        plt.title('Trimer (terminal droplet)',fontsize=60)
    elif(Nclusters==3 and clusterid==1):
        plt.title('Trimer (middle droplet)',fontsize=60)
    plt.ylabel(r'$\frac{{N_{unbonded}}^{t=\infty}}{N_{b}}$',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    fig.tight_layout()
    plt.show()
    plt.close()
    """
    
    #Phase diagram

    fig,ax=plt.subplots(figsize=(21,15),dpi=100)

    x=results[:,0]
    y=results[:,1]
    z=results[:,2]

    x=np.unique(x)
    y=np.unique(y)

    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    X,Y = np.meshgrid(x,y)
    Z=z.reshape(len(y),len(x))

    #print(Z)

    im=plt.pcolormesh(X,Y,Z,cmap='viridis')    

    ax.autoscale_view()
    ax.set_xlabel(r'$R_{droplet}$',fontsize=60)
    ax.set_ylabel(r'$N_{b}$',fontsize=60)

    ax.autoscale_view()

    if(Nclusters==2 and clusterid==0):
        plt.title('Dimer',fontsize=60)
    elif(Nclusters==2 and clusterid==1):
        plt.title('Dimer',fontsize=60)
    elif(Nclusters==3 and clusterid==0):
        plt.title('Trimer (terminal droplet)',fontsize=60)
    elif(Nclusters==3 and clusterid==1):
        plt.title('Trimer (middle droplet)',fontsize=60)

    cb=fig.colorbar(im, ax=ax)
    im.set_clim(0.0,1.0)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(5)
    cb.ax.tick_params(labelsize=50)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    cb.set_label(r'$\frac{{N_{unbonded}}^{t=\infty}}{N_{b}}$',labelpad=+15,fontsize=100)
    fig.tight_layout()
    plt.savefig('final_figures/Nclusters'+str(Nclusters)+'_clusterid'+str(clusterid)+'_RandNpheatmap.png',bbox_inches='tight')
    plt.close()
    
for i in range(0,2,1):
    clusterid=i
    get_phasediagrams(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid)
