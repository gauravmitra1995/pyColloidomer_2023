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

parser=argparse.ArgumentParser()

#Parameters to be passed from outside

#parser.add_argument("--R",type=float)
parser.add_argument("--Np",type=int,default=100)
#parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int,default=3)

#Parameters not to be passed from outside
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


dt=0.001
dimension=2

def round_to_int(val):

    true=val
    nint=int(val)
    if(np.abs(nint-true)<0.5):
        rounded=np.floor(true)
    else:
        rounded=np.ceil(true)
    return np.array(rounded)

def get_phasediagrams(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch):

    epslist=[]
    Rlist=[]
    meandifference_patches_terminaldroplets_dict={}
    results=[]

    for filename in sorted(glob.glob(str(os.getcwd())+"/differencepatches_vs_time_data_new/trimer_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R*_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps*_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".averageallseeds.allruns.difference_terminaldroplets.data"),key=lambda x:(float(os.path.basename(x).split("_")[7].replace('eps','')),float(os.path.basename(x).split("_")[4].replace('R','')))):
        data=np.load(filename,allow_pickle=True)
        eps=float(os.path.basename(filename).split("_")[7].replace('eps',''))
        if(eps not in epslist):
            epslist.append(eps)
        R=float(os.path.basename(filename).split("_")[4].replace('R',''))
        if(R not in Rlist):
            Rlist.append(R)
        meandifferencelist=data[0]
        stddevdifferencelist=data[1]
        label=(R,eps)
        meandifference_patches_terminaldroplets_dict[label]=meandifferencelist
        results.append(np.array([R,eps,meandifferencelist[-1]]))

    results=np.array(results)

    fig,ax=plt.subplots(figsize=(25,18),dpi=100)

    x=results[:,0]
    y=results[:,1]
    z=results[:,2]
    x=np.unique(x)
    y=np.unique(y)
    
    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    X,Y = np.meshgrid(x,y)
    Z=z.reshape(len(y),len(x))

    im=plt.pcolormesh(X,Y,Z,cmap='viridis',shading='nearest')
    
    ax.autoscale_view()
    ax.set_ylabel(r'$\varepsilon$',fontsize=60)
    ax.set_xlabel(r'$R_{\mathrm{droplet}}$',fontsize=60)
    ax.set_title('Trimer: '+r'$N_{b}=$'+str(Np),fontsize=60)

    cb=fig.colorbar(im, ax=ax)
    im.set_clim(0.0,np.max(Z))
    cb.set_label(r'$\frac{|N_{\mathrm{ub},t=t_f}^{(0)}-N_{\mathrm{ub},t=t_f}^{(2)}|}{\frac{1}{2}(N_{\mathrm{ub},t=t_f}^{(0)}+N_{\mathrm{ub},t=t_f}^{(2)})}$',labelpad=+25,fontsize=100)
    cb.ax.tick_params(labelsize=50)
    plt.xticks(np.arange(20,220,20),fontsize=50)
    plt.yticks(np.arange(5.0,30.0,5.0),fontsize=50)
    fig.tight_layout()
    plt.savefig('final_figures/Nclusters3_varyRandeps_asymmetry.png',bbox_inches='tight')
    plt.close()
        
get_phasediagrams(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch)
