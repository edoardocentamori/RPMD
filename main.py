import numpy as np
from modules.rattle import *
from setting import *
from modules.verlet import verlet_algorithm, initialize
from modules.graph import *
from modules.estimators import energy, angular_momentum
import pickle
from matplotlib import pyplot as plt


# file naming 

prepath='/Users/edoardo/Desktop/simulazione_prova/record/'
ID = str(P)+'-'+str(N) + '-' + str(dt) +'-'+ str(id) +'.txt'
path=prepath+ID

q_0,v_0 = initialize()

# Constrain setting 

#fix=[(0,0),(3,0)]
#fixc=[0,1]
#lint=np.sqrt(((q_0[0,0]-q_0[1,0])**2).sum(-1))
#lenf=[((0,0),(1,0),lint)]
#lenfc=[(0,1,lintc)]
#fixc =[0]
fix=None
fixc=None
lenf=None
lenfc=None

u_0=q_0.sum(0)
#lintc=np.sqrt(((u_0[0]-u_0[1])**2).sum(-1))/N


# Actual computation

Q_n,V_n= verlet_algorithm(q_0,v_0,fix=fix,fixc=fixc,lenf=lenf,lenfc=lenfc,norm=1,therm=0)
#Q_n1,V_n1 = verlet_algorithm(q_0,v_0,fix=fix,fixc=fixc,lenf=lenf,lenfc=lenfc,norm=0)

G_n=(Q_n.swapaxes(0,1)).swapaxes(1,2)
#l_n = np.sqrt(((G_n[0,0]-G_n[1,0])**2).sum(-1))

E_n = energy(Q_n,V_n)
#E_n1 = energy(Q_n1,V_n1)
L_n = angular_momentum(Q_n,V_n)

# Data storage

out_file=open(path,'wb+')
pickle.dump([Q_n,V_n,E_n,L_n,dt],out_file)
out_file.close()

# Plotting energy and angular momentum

fig, (ax1,ax2) = plt.subplots(2,1)
x = np.linspace(0,steps*dt,steps-3)
E_n=E_n[3:]
#E_n1=E_n1[3:]
L_n=L_n[3:]

#L_n=l_n[3:]-l_n.mean()
ax1.scatter(x,E_n,c='yellow',s=0.1)
#ax1.scatter(x,E_n1,c='blue',s=0.1)
ax2.scatter(x,L_n,c='red',s=0.1)
ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax2.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.show()