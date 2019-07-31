import numpy as np
import pickle

#q_0=l*np.asarray([[np.cos(2*np.pi*i/N),np.sin(2*np.pi*i/N)] for i in range(N)])

steps = 2000
N = 7
dim = 2

q = np.zeros((N,steps,dim))
u = np.zeros((N,steps,dim))*0j

#C = np.zeros((N,N))

C = np.asarray([[np.exp(-2j*np.pi*i*j/N) for i in range(N)]for j in range(N)])

def omega(i):
    assert i<N, 'index of the mode overflow'
    return 4*np.sin(np.pi*i/N)

def mode(i,a,b):
    sing = np.sin(omega(i)*np.linspace(0,2*np.pi*6, steps))
    cosg = np.cos(omega(i)*np.linspace(0,2*np.pi*6, steps))
    return np.asarray([a*l*(cosg-1j*sing),1j*b*l*(cosg-sing*1j)]).T    

#sin = np.sin(np.linspace(0,2*np.pi*6, steps))
#cos = np.cos(np.linspace(0,2*np.pi*6, steps))
#FREQUENCY MUST BE DIFFERENT!!!

# a^2+b^2=~1
a=1.
b=1.
l=3.
modes = [2,3]

#u[2]=np.asarray([a*l*sin,b*l*sin]).T
#u[3]=np.asarray([b*l*sin,a*l*sin]).T

#u[4]=mode(4,0.2,0.9)
#u[5]=mode(5,0.3,-0.8)
#u[2]=mode(2,a,b)
#u[3]=mode(3,b,a)
#u[6]=mode(6,-0.5,0.5)
u[1]=mode(1,1.,1.)

'''
for m in modes :
    u[m]=np.asarray([a*l*sin,b*l*sin]).T
'''

q= 1/np.sqrt(N)*np.dot(C, u.swapaxes(0,1)).swapaxes(0,1)
q=q.imag
a,b,c=q.shape
q=q.reshape(a,1,b,c)


out_file=open('/Users/edoardo/Desktop/simulazione_prova/pickle_norm.txt','wb')
pickle.dump([q,q,q,q,0],out_file)
out_file.close()
