import numpy as np
from setting import *
import matplotlib.pyplot as plt
import pickle
from modules.estimatorsM import *

# Data extraction

location='/Users/edoardo/Desktop/simulazione_prova/record/'
pathcoord='1-40-0.05-21'+'.txt'
path=location+pathcoord

in_file=open(path, 'rb')
A=pickle.load(in_file)
in_file.close() #aggiunte recentemente, sembrava non necessario

Q_n, V_n, E_n, L_n ,dt= A[0], A[1], A[2], A[3], A[4]

#plotting debug estimators

K = kineticM(V_n)
Pot1 = pot1M(Q_n)
k = time_avarage(K)
p1 = time_avarage(Pot1)
Pot2 = pot2M(Q_n)
p2 = time_avarage(Pot2)


x=np.linspace(0,len(k)*dt,len(k))
fig, ax = plt.subplots()
ax.scatter(x,k,c='red',s=0.1,label='kinetic')
ax.scatter(x,p1,c='blue',s=0.1,label='pot1')
ax.scatter(x,p2,c='green',s=0.1, label='pot2')
ax.legend()
plt.show()