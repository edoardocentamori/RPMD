import numpy as np
import random as rnd
from matplotlib import pyplot as plt
import time 
from matplotlib import animation
import pickle

N = 7
omega2=1.
dt=0.01
steps = 2000
dim = 2
l=3
dt_res=1.

q_0=l*np.asarray([[np.cos(2*np.pi*i/N),np.sin(2*np.pi*i/N)] for i in range(N)])
q_n=q_0+(np.random.rand(N,dim)*0.8-0.4)

v_n=np.random.rand(N,dim)*0.2-0.1

def dV(q):
    V=[2*q[0]-q[1]-q[-1]]
    for i in range(1,N-1):
        V.append(2*q[i]-q[i+1]-q[i-1])
    V.append(2*q[-1]-q[0]-q[-2])
    return omega2*np.asarray(V)

def verlet_step(q,v):
    q1=q+v*dt
    lamda_v=2/dt*v[0]-dV(q)[0]
    #lamda = 2*(1.-omega2*dt**2)*q1[0]-q[0]+dt**2*omega2*(q1[1]+q1[-1])
    lamda_r=(q1[0]/dt+v[0]-dt*dV(q1)[0])/(dt/2)-lamda_v
    v12=v-dt/2*dV(q)
    v12[0]-=dt/2*lamda_r
    #v12[0]-=lamda*dt/2
    v1=v12-dt/2*dV(q)
    v1[0]-=dt/2*lamda_v
    return q1,v1
    
q_nm1, v_nm1 = q_n, v_n
q_n, v_n = verlet_step(q_n,v_n)
    
def verlet_algorithm(q_0,v_0):
    Q_n, V_n= [], []
    q_n, v_n = q_0, v_0
    for _ in range(steps):
        Q_n.append(q_n)
        V_n.append(v_n)
        q_n, v_n = verlet_step(q_n, v_n)
    return np.asarray(Q_n), np.asarray(V_n)
    

Q_n, V_n = verlet_algorithm(q_n,v_n)
Q_nc, V_nc = np.asarray([Q_n.sum(1)])/N, np.asarray([V_n.sum(1)])/N
Q_nc, V_nc = Q_nc.swapaxes(0,1), V_nc.swapaxes(0,1)
Q_n = np.concatenate((Q_n,Q_nc), 1)
V_n = np.concatenate((V_n,V_nc), 1)


out_file=open('/Users/edoardo/Desktop/simulazione_prova/pickle_sym.txt','wb')
pickle.dump([Q_n,V_n],out_file)
out_file.close()


