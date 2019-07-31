import numpy as np
import time
from setting import *
from modules.verlet import initialize
#start = time.time()


#N = 20000
#(N,P,dim)

assert N%2==0, 'N must be even'

q_0, v_0 = initialize()

def C_matrix():
    '''
    j-k inverted w.r.t the definition in the paper (C = C_kj)
    took 2.3 s for N=10000
    '''
    C = np.zeros((N,N))
    j=np.arange(N)
    one= np.ones(N)
    for k in range(N):
        if k==0:
            C[k]=np.sqrt(1/N)*one
        elif 1<=k<=N/2-1:
            C[k]=np.sqrt(2/N)*np.cos(2*np.pi*j*k/N)
        elif k==N/2:
            C[k]=np.sqrt(1/N)*(-1)**(j)
        elif N/2+1<=k<=N-1:
            C[k]=np.sqrt(2/N)*np.sin(2*np.pi*j*k/N)
    return C

C=C_matrix()
q_1=np.dot(C,q_0.swapaxes(0,1)) #works, now you have to get correct the multiplication stuff
q_2=np.dot(C.T,q_1.swapaxes(0,1))
#end = time.time()
#print(end - start)

'''
for j in range(N):
    for k in range(N):
        if k==0:
            C[j,k]=np.sqrt(1/N)
        elif 1<=k<=N/2-1:
            C[j,k]=np.sqrt(2/N)*np.cos(2*np.pi*j*k/N)
        elif k==N/2:
            C[j,k]=np.sqrt(1/N)*(-1)**(j)
        elif N/2+1<=k<=N-1:
            C[j,k]=np.sqrt(2/N)*np.sin(2*np.pi*j*k/N)
#this one took 310 s for N=10000
'''