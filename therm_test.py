import numpy as np
from setting import *
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
from modules.normal import v_norm
from modules.estimators import mode_energy
import pickle
from scipy.optimize import curve_fit

location='/Users/edoardo/Desktop/simulazione_prova/record/'
pathcoord='1-10-0.005-16'+'.txt'
path=location+pathcoord

in_file=open(path, 'rb')
A=pickle.load(in_file)
in_file.close() #aggiunte recentemente, sembrava non necessario
Q_n, V_n, E_n, L_n ,dt= A[0], A[1], A[2], A[3], A[4]


U_n = []
VU_n = []
for q in Q_n:
    U_n.append(v_norm(q))
    
for v in V_n:
    VU_n.append(v_norm(v))
    
U_n = np.asarray(U_n)
VU_n = np.asarray(VU_n)

bins=20

#a = np.random.normal(size=(10000))

E = mode_energy(U_n,VU_n).swapaxes(0,1)


'''
# here is debugging
fig5, ax5 = plt.subplots()
for i in range(N):
    E1=E[i]
    x=np.arange(len(E1))
    ax5.scatter(x,E1,s=0.01)
#E1=E[10]
#x=np.arange(len(E1))
#ax5.scatter(x,E1,s=0.01)
ax5.set_title('energy')
#ax5.set_xlabel('modes')
#ax5.set_ylabel('fitted beta')

plt.show()

u=(U_n.swapaxes(0,1)[10]).swapaxes(0,-1)[0,0]
T=np.linspace(0,steps*dt,steps)
fig6, ax6 = plt.subplots()
ax6.scatter(T,u,c='b',s=0.1)
plt.show()

#here is end of debugging
'''



E1 = E[5]
#E2 = E[100]
#E3 = E[250]


fig, ax = plt.subplots()

ax.hist(E1,bins,density=True,facecolor='b',alpha=0.75)
#ax.hist(E2,bins,density=True,facecolor='g',alpha=0.75)
#ax.hist(E3,bins,density=True,facecolor='y',alpha=0.75)

#ax.hist(statistic,bin_edges,color='blue')
#ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
#ax.legend()
plt.show()



# the histogram of the data
#n, bins, patches = plt.hist(x, 50, density=True, facecolor='g', alpha=0.75)


plt.xlabel('dummy variable')
plt.ylabel('Probability distribution')
#plt.title('Histogram of IQ')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(40, 160)
#plt.ylim(0, 0.03)
plt.grid(True)
plt.show()


statistic, bin_edges, bin_number = binned_statistic(E1, np.arange(bins), statistic = 'count', bins=bins)



x=np.zeros(len(bin_edges)-1)

for i in range(len(bin_edges)-1):
    x[i]=(bin_edges[i]+bin_edges[i+1])/2
print(x)

y= statistic/len(E1)
#y=y[1:]

#x=x[1:]

#funzione di fit
def F(x,C,b):
    return C*x**(dim-1)*np.exp(-b*x)
    
init=(1.,beta/N)
popt,pcov=curve_fit(F,x,y,init)


print("C = ", popt[0], " +/- ", np.sqrt(pcov[0,0]) )
print("beta = ", popt[1], " +/- ", np.sqrt(pcov[1,1]))

x1=np.linspace(0,np.max(E1),10000)
y1=F(x1,popt[0],popt[1])


fig1, ax1 = plt.subplots()

ax1.scatter(x,y,c='b')
ax1.scatter(x1,y1,c='g',s=0.01)
#ax1.set_yscale('log')


plt.show()

fitted_beta=[]



def find_beta(E):
    statistic, bin_edges, bin_number = binned_statistic(E, np.arange(bins), statistic = 'count', bins=bins)
    x=np.zeros(len(bin_edges)-1)
    for i in range(len(bin_edges)-1):
        x[i]=(bin_edges[i]+bin_edges[i+1])/2
    y= statistic/len(E)
    y=y[1:]
    x=x[1:]
    popt,pcov=curve_fit(F,x,y,init)
    return popt[1]
    

for i in range(N):
    fitted_beta.append(find_beta(E[i]))


fig2, ax2 = plt.subplots()
modes=np.arange(N)

ax2.scatter(modes,fitted_beta,c='b')
ax2.set_xlabel('modes')
ax2.set_ylabel('fitted beta')

plt.show()
