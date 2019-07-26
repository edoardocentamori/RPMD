import numpy as np
from matplotlib import pyplot as plt

a=3000
N=100

def add(i,N,a):
    assert i<N, 'problem in add'
    return 1./(1+(a*np.sin(np.pi*i/N)**2))
    
def sum(N,a):
    cum = 0.
    for i in range(N):
        #np.exp(2*np.pi*1j/N)
        cum +=add(i,N,a)
    print(cum/N)
    return cum

#y=np.asarray([sum(j+1) for j in range(N-1)])

y=np.asarray([sum(10,b) for b in range(a)])
#yre=y.real
#yim=y.imag
x=np.linspace(0,a,a)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x,y)
plt.yscale('log')
plt.xscale('log')
#ax.plot(x,yre)
#ax.plot(x,yim)
plt.show()
