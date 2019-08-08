import numpy as np
import matplotlib.pyplot as plt
from setting import *
from modules.verlet import verlet_algorithm, zero_init
from modules.estimators import *

q_0,v_0 = zero_init()


'''
final = []

for beta in [0.1,0.2,0.4,0.6,0.8,1.,2.,3.,4.]:
    Q_n,V_n = verlet_algorithm(q_0,v_0,beta,fix=fix,fixc=fixc,lenf=lenf,lenfc=lenfc,norm=1,therm=1,debug=False)
    E1 = primitive_energy(Q_n,V_n)
    E2 = real_primitive_energy(Q_n,V_n)
    e1 = time_avarage(E1)
    e2 = time_avarage(E2)
    final.append((e1,e2))
'''

x = np.linspace(0,steps*dt,steps)
fig, ax = plt.subplots()
for e1,e2 in final:
    ax.scatter(x,e1,c='blue',s=0.1,label='primitive')
    ax.scatter(x,e2,c='green',s=0.1,label='real prim')
ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax.set_title('Time avarage')
ax.legend()
plt.show()

