from setting import *
from modules.normal import v_norm, v_stand, w
from numpy.random import normal
import numpy as np
import matplotlib.pyplot as plt

gamma = 2*w
gamma[0] = 1/t0

c1=np.exp(-dt/2*gamma)
c2=np.sqrt(1-c1**2)

#I'm adding the beta -> beta/N in the thermalization, let's hope it save the day

def thermal_step(v,beta):
    '''
    It might have some problem for P>1 must fix
    '''
    vu = v_norm(v)
    #print(vu[5])
    #print(c1[5])
    vu = (c1*vu.swapaxes(0,-1)).swapaxes(0,-1)
    #print(vu[5])
    #print('----')
    # to remove the beta/N just sobstitute N-> 1 in the sqrt
    # removed for the moment
    #print(np.sqrt(N/beta)*c2[5])
    ran = normal(size=(dim,P,N))
    #print(ran.swapaxes(0,-1)[5,0])
    delta = (np.sqrt(N/beta)*c2*ran).swapaxes(0,-1)
    #print(delta[5,0])
    vu+= delta
    #print(vu[5])
    #print('-----')
    return v_stand(vu) 
    
def debug_thermal_step(vu,beta):
    vu = (c1*vu.swapaxes(0,-1)).swapaxes(0,-1)
    ran = normal(size=(dim,P,N))
    delta = (np.sqrt(N/beta)*c2*ran).swapaxes(0,-1)
    vu+= delta
    return vu 

fig, ax = plt.subplots()
x = np.linspace(0,N,N)

ax.scatter(x,c1,c='green',s=0.1,label='c1')
ax.scatter(x,c2,c='blue',s=0.1, label='c2')
ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax.set_xlabel('mode')
ax.set_ylabel('C1')
ax.legend()
ax.set_title('C1 and C2')
plt.show()