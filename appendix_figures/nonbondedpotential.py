import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("poster", font_scale=1.0)

radius=1.0
rcut=2*radius
r=np.linspace(0,2*rcut,1000)
ron=0.1*rcut
epsilon=500.0

#Details about the soft quartic potential plotted below are described in Appendix under "Non-bonded interaction details"

fig=plt.figure(figsize=(8,6))
U=np.zeros(r.size)
S=np.zeros(r.size)
V=np.zeros(r.size)
for i in range(r.size):
    if(r[i]<ron):
        S[i]=1
    elif(ron<=r[i] and r[i]<=rcut):
        S[i]=((rcut**2-r[i]**2)**2)*(rcut**2 + 2*r[i]**2 - 3*ron**2)/((rcut**2-ron**2)**3)
    elif(r[i]>rcut):
        S[i]=0
    if(r[i]<rcut):
        U[i]=epsilon*(1-(r[i]/rcut)**4)
    else:
        U[i]=0
    shift=0
    if(ron<rcut):
        V[i]=S[i]*U[i]
    elif(ron>=rcut):
        V[i]=U[i]-shift
p=plt.plot(r,V,color='red',label=r'$V(r) = U(r)S(r)$',linewidth=5.0)
plt.plot(r,U,color='darkgreen',label=r'$V(r) = U(r)$',linestyle='--',linewidth=5.0)
plt.xlabel(r'$r$',size=30)
plt.ylabel(r'$V(r)$',fontsize=30)
plt.axvline(x=rcut,color='black',linestyle='--',linewidth=3.0)
plt.legend(loc='upper right',fontsize=20)
fig.tight_layout()
plt.savefig('softpotential.png',bbox_inches='tight')
plt.close()
