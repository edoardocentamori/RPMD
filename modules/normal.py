import numpy as np
from setting import *
    
#j-k inverted w.r.t the definition in the paper (C = C_kj)
#took 2.3 s for N=10000

assert N%2==0, 'N must be even'
C = np.zeros((N,N))
j=np.arange(N)
for k in range(N):
    if k==0:
        C[k]=np.sqrt(1/N)*np.ones(N)
    elif 1<=k<=N/2-1:
        C[k]=np.sqrt(2/N)*np.cos(2*np.pi*j*k/N)
    elif k==N/2:
        C[k]=np.sqrt(1/N)*(-1)**(j)
    elif N/2+1<=k<=N-1:
        C[k]=np.sqrt(2/N)*np.sin(2*np.pi*j*k/N)

w=2*np.sqrt(omega2)*np.sin(np.pi/N*np.arange(N))



numerical_sin_w = np.zeros(N)
# solving 0/0 problem
for i in range(N):
    if w[i] != 0:
        numerical_sin_w[i] = 1/w[i]*np.sin(w[i]*dt)
    else :
        numerical_sin_w[i] = dt
        
def get_norm(q,v):
    return np.dot(C,q.swapaxes(0,1)), np.dot(C,v.swapaxes(0,1))

def v_norm(v):
    return np.dot(C,v.swapaxes(0,1))
    
def prop_norm(u,vu):
    vu1 = np.cos(w*dt)*vu.swapaxes(0,-1)-w*np.sin(w*dt)*u.swapaxes(0,-1)
    u1 = np.cos(w*dt)*u.swapaxes(0,-1)+numerical_sin_w*vu.swapaxes(0,-1)
    return u1.swapaxes(0,-1), vu1.swapaxes(0,-1)

def get_stand(u,vu):
    return np.dot(C.T,u.swapaxes(0,1)), np.dot(C.T,vu.swapaxes(0,1))
 
def v_stand(vu):
    return np.dot(C.T,vu.swapaxes(0,1))

    