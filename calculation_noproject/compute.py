from setting import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

a=4*N**2/(beta**2*ext_omega2)

def add(i,N,a):
    assert i<N, 'problem in add'
    return 1./(1+a*(np.sin(np.pi*i/N))**2)
    
def sum(N,a):
    cum = 0.
    for i in range(N):
        cum +=add(i,N,a)
    return cum



expected1 = 1/(beta*ext_omega2)*sum(N,a)/2

#this is the expected value for pot1

print(expected1)

theorical_expected1 = 1/2*np.sqrt(ext_omega2)*beta/2*1/np.tanh(beta*np.sqrt(ext_omega2)/2)

print(theorical_expected1)


'''
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x,y)
plt.yscale('log')
plt.xscale('log')
#ax.plot(x,yre)
#ax.plot(x,yim)
plt.show()
'''

'''
N=100
beta=1

energies =[]

WE = np.arange(0.01,0.3,0.001)

c=4*N**2/beta**2 #added N -> N**2 to test 
print(c)

for we in WE:
    a =c/we**2
    energies.append(sum(N,a))
    


#funzione di fit
def F(x,a,b):
    return a*x/2/np.tanh(b*x/2)
    
init=(1.,1.)
popt,pcov=curve_fit(F,WE,energies,init)


print("a = ", popt[0], " +/- ", np.sqrt(pcov[0,0]) )
print("b = ", popt[1], " +/- ", np.sqrt(pcov[1,1]))

x1=np.linspace(0.01,np.max(WE),10000)
y1=F(x1,popt[0],popt[1])


fig, ax = plt.subplots()
ax.scatter(WE,energies,c='red',s=0.1)
ax.scatter(x1,y1,c='green',s=0.01)
plt.show()

#NOW with N-> N**2 it works, and fit the correct thing, the problem is that w**2 = N /beta**2 and not N**2/beta**2
    
'''