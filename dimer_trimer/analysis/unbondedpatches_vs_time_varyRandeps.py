from functions_cluster_analysis import *
from scipy.optimize import curve_fit

parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside
#parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int)
#parser.add_argument("--R",type=float)
parser.add_argument("--Np",type=int,default=100)

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


def fit_double_exponential(meanfractionlist,timelist):
    a=meanfractionlist[100]
    b=meanfractionlist[-1]
    k1=1.0/timelist[100]
    k2=1.0/timelist[-1]

    fit_params,pcov=curve_fit(lambda t,a,k1,b,k2:  (meanfractionlist[0]-a)*np.exp(-k1*timelist)+(a-b)*np.exp(-k2*timelist)+b, timelist, meanfractionlist, p0=(a,k1,b,k2),maxfev=10000, bounds=(0,np.inf))

    a=fit_params[0]
    k1=fit_params[1]
    b=fit_params[2]
    k2=fit_params[3]
    sat=b
    truesat=meanfractionlist[-1]
    
    recruitment_time=1.0/np.max((k1,k2))
    patchsaturationtime=1.0/np.min((k1,k2))
    function=(meanfractionlist[0]-a)*np.exp(-k1*np.array(timelist))+ (a-b)*np.exp(-k2*np.array(timelist)) +b

    return sat,truesat,recruitment_time,patchsaturationtime,function


def generate_plots_and_data(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid):

    timesteplist=[]
    labellist=[]
    epslist=[]
    Rlist=[]
    meanunbondedbindersdict={}
    stddevunbondedbindersdict={}


    #for dimer, terminal 1 and terminal 2 droplets should essentially have identical data 
    if(Nclusters==2):
        simulationtype='dimer'
        if(clusterid==0):
            droplet='terminal 1'
        elif(clusterid==1):
            droplet='terminal 2'
    elif(Nclusters==3):
        simulationtype='trimer'
        if(clusterid==0):
            droplet='terminal 1'
        elif(clusterid==1):
            droplet='middle'
        elif(clusterid==2):
            droplet='terminal 2'

    for filename in sorted(glob.glob(str(os.getcwd())+"/unbondedpatches_vs_time_data/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R*_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps*_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".averageallseeds.allruns.num_unbondedpatches.data"),key=lambda x:(float(os.path.basename(x).split("_")[4].replace('R','')),float(os.path.basename(x).split("_")[7].replace('eps','')))):
        data=np.load(filename,allow_pickle=True)
        R = float(os.path.basename(filename).split("_")[4].replace('R',''))
        eps = float(os.path.basename(filename).split("_")[7].replace('eps',''))
        label=(R,eps)

        meanunbondedbindersdict[label]=data[0]
        stddevunbondedbindersdict[label]=data[1]
        l=len(data[0])
        labellist.append(label)
        if(eps not in epslist):
            epslist.append(eps)
        if(R not in Rlist):
            Rlist.append(R)
        timesteplist=data[2]

    timelist=np.zeros(len(timesteplist))
 
    for i in range(len(epslist)):
        number_of_plots=len(Rlist)
        eps=epslist[i]
        fig,ax=plt.subplots(figsize=(20,15),dpi=100)
        #colormap=plt.cm.jet
        #colors = [colormap(i) for i in np.linspace(0, 1,number_of_plots)]
        colors=['blue','darkorange','green','crimson','darkmagenta','teal','deeppink','limegreen']
        #ax.set_prop_cycle('color', colors)

        saturationfractionlist=[] 
        recruitment_time_list=[]
        patchsaturation_time_list=[]

        r=0
        for j in range(len(Rlist)):
            R=Rlist[j]
            label=(R,eps)       
            meanfractionlist=np.array(np.array(meanunbondedbindersdict[label][clusterid])/(Np))
            stddevfractionlist=np.array(np.array(stddevunbondedbindersdict[label][clusterid])/(Np))
           
            timesteplist=np.array(timesteplist)
            timelist=timesteplist*dt
            p=ax.plot(timelist[::5],meanfractionlist[::5],marker='.',markersize=18,color=colors[r],label=r'$R = $'+str(label[0]))

            #calculate the saturation fraction of unbonded binders according to the regime of epsilon (do double exponential fitting for high epsilon) 
            if(eps>=13.8):
                sat,truesat,recruitment_time,patchsaturationtime,function=fit_double_exponential(meanfractionlist,timelist)
                if(sat>=0.0 and sat<=1.0 and np.abs(sat-truesat)<=0.1):
                    saturationfractionlist.append(np.round(sat,4))
                else:
                    saturationfractionlist.append(np.round(truesat,4))
                recruitment_time_list.append(recruitment_time)
                patchsaturation_time_list.append(patchsaturationtime)
            else:
                sat=np.mean(meanfractionlist[1000:],axis=0)
                saturationfractionlist.append(np.round(sat,4))


            #plot the fitting functions according to the regime of epsilon 
            if(eps>=13.8):
                ax.plot(timelist,function,color=p[0].get_color(),linewidth=5.0)
            else:
                ax.plot(timelist[10:],np.ones(timelist[10:].shape[0])*sat,color=p[0].get_color(),linewidth=5.0)
            r+=1

        print("epsilon:",eps," droplet:",droplet,saturationfractionlist)

        plt.ylabel(r'$\frac{N_{unbonded}}{N_{b}}$',labelpad=20,fontsize=100)
        plt.xlabel("Simulation time (in HOOMD units)",labelpad=20,fontsize=60)
        
        if(Nclusters==2 and clusterid==0):
            plt.title('Dimer: '+r'$\varepsilon = $'+str(eps),fontsize=60)
        elif(Nclusters==2 and clusterid==1):
            plt.title('Dimer: '+r'$\varepsilon = $'+str(eps),fontsize=60)
        elif(Nclusters==3 and clusterid==0):
            plt.title('Trimer (terminal droplet): '+r'$\varepsilon = $'+str(eps),fontsize=60)
        elif(Nclusters==3 and clusterid==1):
            plt.title('Trimer (central droplet): '+r'$\varepsilon = $'+str(eps),fontsize=60)
        
        plt.xticks(np.arange(0,210000,50000),fontsize=50)
        plt.yticks(fontsize=50)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(5)
        plt.legend(loc='upper right',ncol=2,prop={'size': 40})
        fig.tight_layout()
        plt.savefig('final_figures/unbondedbinders_vs_time'+'_Nclusters'+str(Nclusters)+'_clusterid'+str(clusterid)+'_Np'+str(Np)+'_epsilon'+str(eps)+'_varyR.png',bbox_inches='tight')
        plt.close()
        
        #for each epsilon, save the saturation fraction list corresponding to each of the droplet radii
        output_saturationfraction=os.path.join(os.getcwd(),"saturation_data","saturationfraction_restl"+str(r0)+"_Nc"+str(Nclusters)+"_clusterid"+str(clusterid)+"_Np"+str(Np)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_eps"+str(eps)+"_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".varyingRandeps.data")
        x=np.array([Rlist,saturationfractionlist],dtype=object)
        x.dump(output_saturationfraction)
            
for i in range(0,2,1):
    clusterid=i
    generate_plots_and_data(Nclusters,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid)



