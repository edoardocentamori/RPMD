import numpy as np
from modules.rattle import *
from setting import *
from modules.verlet import verlet_algorithm, initialize, zero_init
from modules.verletM import verlet_algorithmM
from modules.graph import *
from modules.estimators import *
from modules.estimatorsM import primitive_energyM, real_primitive_energyM,real_primitive_energyH2M, KH2M, pot1H2M,pot2H2M
import pickle
from matplotlib import pyplot as plt
import time
from modules.estimatorsM import energyM, angular_momentumM
from modules.normal import eta

# Main is used to run the entire model and sketch some basic graph used mainly for debugging purpouses, output data are saved as byte files and processed elsewere

# This part defines the output file in a format that incorporate the data contained in it

prepath='/Users/edoardo/Desktop/simulazione_prova/record/'
ID = str(P)+'-'+str(N) + '-' + str(dt) +'-'+ str(id) +'.txt'
path=prepath+ID


#q_0,v_0 = initialize()
q_0,v_0 = zero_init()


# Constrain setting, some random constrain to check functionality
# Depending on what you need to test you can use different ones [debugging purpouse]


#fix=[(0,0),(3,0)]
#fixc=[0]
#lint=np.sqrt(((q_0[0,0]-q_0[1,0])**2).sum(-1))
#lenf=[((0,0),(1,0),lint)]
#u_0=q_0.sum(0)
#lintc=np.sqrt(((u_0[0]-u_0[1])**2).sum(-1))/N
#lenfc=[(0,1,lintc)]
#fixc =[0]
#fixcm = [0,1]

# These flags tells you which constraint are enabled

fixL=True
fix=None
fixc=None
fixcm=None
lenf=None
lenfc=None


#mtot= m.sum()
#linear = M
#linear = np.zeros((N,P,dim))
#linear[:,0,:]=np.ones((N,dim))
linear = None

# Actual computation, this function run the entire program, the time of the operation is monitored to check for efficiency

start = time.time()

Q_n,V_n= verlet_algorithmM(q_0,v_0,beta,fix=fix,fixc=fixc,fixcm=fixcm,lenf=lenf,lenfc=lenfc,norm=1,therm=True,debug=False,H2=True,linear=linear,fixL=fixL)

end = time.time()
print(end - start)

#Here are defined several observables of interest in the system

#G_n=(Q_n.swapaxes(0,1)).swapaxes(1,2)
#l_n = np.sqrt(((G_n[0,0]-G_n[1,0])**2).sum(-1))

E_n = energyM(Q_n,V_n)
#L_n = angular_momentumM(Q_n,V_n)
L_n= np.zeros(E_n.size) #d1 dimension

# new estimators

#E1 = fictitious_energy(Q_n,V_n)
E2 = primitive_energyM(Q_n,V_n)
E3 = virial_energy(Q_n,V_n)
E4 = real_primitive_energyM(Q_n,V_n)
E5 = real_primitive_energyH2M(Q_n,V_n)
k = KH2M(Q_n,V_n)
p1 = pot1H2M(Q_n,V_n)
p2 = pot2H2M(Q_n,V_n)


# Data storage

out_file=open(path,'wb+')
pickle.dump([Q_n,V_n,E_n,L_n,dt],out_file)
out_file.close()

#From here on there are only a lot of graphs used mainly for debugging and informative purpouses.

# Plotting energy and angular momentum

fig, (ax1,ax2) = plt.subplots(2,1)
x = np.linspace(0,steps*dt,steps-3)
E_n=E_n[3:]
L_n=L_n[3:]

#L_n=l_n[3:]-l_n.mean()
ax1.set_title('Energy')
ax1.scatter(x,E_n,c='yellow',s=0.1)
ax1.set_ylabel('Energy [a.u.]')
ax2.scatter(x,L_n,c='red',s=0.1)
ax2.set_title('Angular momentum')
ax2.set_ylabel('Angular momentum[a.u.]')
ax2.set_xlabel('steps [thousands]')
ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax2.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.show()


# plotting energy estimators

fig1, ax3 = plt.subplots()
#E1=E1[3:]
E2=E2[3:]
E3=E3[3:]
E4=E4[3:]
E5=E5[3:]
if P==2:
    k = k[3:]
    p1 = p1[3:]
    p2 = p2[3:]


#ax3.scatter(x,E1,c='yellow',s=0.1,label='fictitious')
ax3.scatter(x,E2,c='blue',s=0.1,label='primitive')
ax3.scatter(x,E3,c='red',s=0.1,label='virial')
ax3.scatter(x,E4,c='green',s=0.1,label='real prim')
ax3.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax3.legend()
plt.show()


#plotting avaraged estimators

fig2, ax4 = plt.subplots()

e2=time_avarage(E2)
e3=time_avarage(E3)
e4=time_avarage(E4)
e5=time_avarage(E5)
if P==2:
    k=time_avarage(k)
    p1=time_avarage(p1)
    p2=time_avarage(p2)

x1 = np.linspace(0,steps,steps-3)

ax4.scatter(x1,e2,c='blue',s=0.1,label='primitive')
#ax4.scatter(x,e3,c='red',s=0.1,label='virial')
ax4.scatter(x1,e4,c='green',s=0.1,label='real prim')
ax4.set_ylabel('Energy [a.u.]')
ax4.set_xlabel('steps')



ax4.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax4.set_title('Time avarage')
ax4.legend()
plt.show()

#fE4 = E4[20000:]
#fe4 = time_avarage(fE4)
#print(fe4[-1]) 

#check conservation of constrain

u = (Q_n.swapaxes(0,2)[0].sum(0)/N).swapaxes(0,1)[0]
x=np.linspace(0,steps*dt,steps)
fig3, ax5= plt.subplots()

ax5.scatter(x,u,c='blue',s=0.1)
ax5.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax5.set_title('Centroid position')
ax5.set_ylabel('Centroid position [a.u.]')
ax5.set_xlabel('steps')
plt.show()

u1 = (Q_n.swapaxes(0,2)[0].sum(0)/N).swapaxes(0,1)
u2 = (Q_n.swapaxes(0,2)[1].sum(0)/N).swapaxes(0,1)

d = ((u1-u2)**2).sum(0)

fig4, ax6= plt.subplots()

#ax6.set_ylim(221.129,221.147)
ax6.scatter(x,d,c='green',s=0.1)
ax6.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax6.set_title('Length between the centroids')
ax6.set_ylabel('Length between the centroids [a.u.]')
ax6.set_xlabel('steps [thousands]')

plt.show()

C_n = Q_n.sum(1)/N
Cm_n = (m[0]*C_n[:,0,:]+m[1]*C_n[:,1,:])/(m[0]+m[1])

Cm_n = Cm_n[:,0] #chosen axes x to show

x=np.linspace(0,steps*dt,steps)*20
fig5, ax7= plt.subplots()

ax7.scatter(x,Cm_n,c='blue',s=0.1)
ax7.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax7.set_title('Center of mass position')
ax7.set_ylabel('C.M. position [a.u.]')
ax7.set_xlabel('steps')
plt.show()


fig6, ax8= plt.subplots()
x1 = np.linspace(0,steps,steps-3)
ax8.scatter(x1,e5,c='blue',s=0.1,label='primitive H2')
#ax8.scatter(x1,k,c='red',s=0.1,label='kinetic')
#ax8.scatter(x1,p1,c='blue',s=0.1,label='p1')
#ax8.scatter(x1,p2,c='green',s=0.1,label='p2')
ax8.set_title('Primitive energy H2')
ax8.set_ylabel('primitive energy [a.u.]')
ax8.set_xlabel('steps')
ax8.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax8.legend()
plt.show()



C_n = Q_n.sum(1)/N
meta = C_n[:,0,:]-C_n[:,1,:]
modmeta = np.sqrt((meta**2).sum(1))
print(modmeta)
for i in range(dim):
    meta[:,i] = meta[:,i]/modmeta
print(meta)
#meta /= np.sqrt((meta**2).sum(1))
dif= ((meta-eta)**2).sum(1)


x=np.linspace(0,steps*dt,steps)

fig7, ax9= plt.subplots()

ax9.scatter(x,dif,c='blue',s=0.1,label='deviation')
ax9.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax9.set_title('deviation from the right direction')
ax9.set_ylabel('Deviation [a.u.]')
ax9.set_xlabel('steps')
plt.show()
#the discrete oscillations are due to C1 and C2 swhiching place in a oscillation
