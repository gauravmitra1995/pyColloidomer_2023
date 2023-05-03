from functions_cluster_analysis import *
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression


parser=argparse.ArgumentParser()

#Parameters REQUIRED to be passed as arguments from outside
parser.add_argument("--epsilon")
parser.add_argument("--Nclusters",type=int,default=2)
parser.add_argument("--R",type=float,default=50.0)
parser.add_argument("--Np",type=int,default=100)
parser.add_argument("--gammapatch",type=float)


#Parameters not mandatory to be passed as arguments from outside
parser.add_argument("--metropolis",default=1,type=int)
parser.add_argument("--kAB",type=float,default=200.0)
parser.add_argument("--radiusB",type=float,default=1.0)
parser.add_argument("--kspring",type=float,default=10.0)
parser.add_argument("--dimension",type=int,default=2)
parser.add_argument("--gammaA",type=float,default=0.1)
parser.add_argument("--r0",type=float,default=2.0)
parser.add_argument("--dt",default=0.001,type=float)

args = parser.parse_args()
locals().update(vars(args))

if(epsilon!='infinite'):
    epsilon=float(epsilon)
else:
    epsilon=str(epsilon)


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


def fit_single_exponential(meanfractionlist,timelist):
    b=meanfractionlist[-1]
    k=1.0/timelist[-1]

    fit_params,pcov=curve_fit(lambda t,b,k:  (meanfractionlist[0]-b)*np.exp(-k*timelist)+b, timelist, meanfractionlist, p0=(k,b),maxfev=10000, bounds=(0,np.inf))

    b=fit_params[0]
    k=fit_params[1]
    sat=b
    truesat=meanfractionlist[-1]

    recruitment_time=1.0/(k)
    function= (meanfractionlist[0]-b)*np.exp(-k*np.array(timelist)) + b 

    return sat,truesat,recruitment_time,function



def generate_plots_and_data(Nclusters,Np,R,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid):

    timesteplist=[]
    labellist=[]
    konlist=[]
    kofflist=[]

    #for dimer, terminal 1 and terminal 2 droplets should essentially have identical data
    if(Nclusters==2):
        simulationtype='dimer'
        if(clusterid==0):
            droplet='terminal 1'
        elif(clusterid==1):
            droplet='terminal 2'
 
    saturationfractionlist=[]
    recruitment_time_list=[]
    patchsaturation_time_list=[]

    recruitment_time_dict={}

    for filename in sorted(glob.glob(str(os.getcwd())+"/data_allseeds/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/epsilon"+str(epsilon)+"/"+str(simulationtype)+"_restl"+str(r0)+"_Nc"+str(Nclusters)+"_Np"+str(Np)+"_R"+str(R)+"_rB"+str(radiusB)+"_kAB"+str(kAB)+"_kon*_koff*_dim"+str(dimension)+"_kspring"+str(kspring)+"_gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+".allsamples.allruns.num_unbondedpatches.data"),key=lambda x:(float(os.path.basename(x).split("_")[7].replace('kon','')),float(os.path.basename(x).split("_")[8].replace('koff','')))):
        data=np.load(filename,allow_pickle=True)
        kon = float(os.path.basename(filename).split("_")[7].replace('kon',''))
        koff = float(os.path.basename(filename).split("_")[8].replace('koff',''))
        label=(kon,koff)
        unbondedbinders=data[0]
        labellist.append(label)

        if(kon not in konlist):
            konlist.append(kon)
        if(koff not in kofflist):
            kofflist.append(koff)

        timesteplist=data[1]
        timelist=np.zeros(timesteplist.shape[0])
        timesteplist=np.array(timesteplist)
        timelist=timesteplist*dt
        fractions=np.zeros(unbondedbinders.shape)

        saturationfractionlist=[]
        recruitment_time_list=[]
        patchsaturation_time_list=[]

        for i in range(unbondedbinders.shape[0]):
            ub=unbondedbinders[i]
            fractions[i]=ub/Np
            if(epsilon>13):
                sat,truesat,recruitment_time,patchsaturationtime,function=fit_double_exponential(fractions[i],timelist)
                if(sat>=0.0 and sat<=1.0 and np.abs(sat-truesat)<=0.1):
                    saturationfractionlist.append(np.round(sat,4))
                else:
                    saturationfractionlist.append(np.round(truesat,4))
                recruitment_time_list.append(recruitment_time)
                patchsaturation_time_list.append(patchsaturationtime)
            else:
                sat,truesat,recruitment_time,function=fit_single_exponential(fractions[i],timelist)
                if(sat>=0.0 and sat<=1.0 and np.abs(sat-truesat)<=0.1):
                    saturationfractionlist.append(np.round(sat,4))
                else:
                    saturationfractionlist.append(np.round(truesat,4))
                recruitment_time_list.append(recruitment_time)

        recruitment_time_dict[label]=recruitment_time_list

    output_file_tau=str(os.getcwd())+"/data_allseeds/gammaA"+str(gammaA)+"_gammapatch"+str(gammapatch)+"/epsilon"+str(epsilon)+"/dimer_epsilon"+str(epsilon)+"_gammapatch"+str(gammapatch)+"_R"+str(R)+".allsamples.recruitmenttimes.txt"
    x=recruitment_time_dict
    with open(output_file_tau, "w") as file:
        for key, value in x.items():
            file.write(f"{key}: {value}\n")
    
for i in range(0,1,1):
    clusterid=i
    generate_plots_and_data(Nclusters,Np,R,r0,kAB,radiusB,dimension,kspring,gammaA,gammapatch,clusterid)



