import numpy as np
from matplotlib import pyplot as plt
import pickle

P = 2 # masses
N = 5 # beads per mass

steps = 20000
dim = 2
# time scales
dt=0.001
dt_res=1. 
# lenght scales
large_scale=20. # scale in which masses are randomly sampled
disp_scale=0.8 # scale in which beads get moved from lattice randomly
l=3 # initial distance beads from centroid
# velocity scales
vel_beads_scale = 0.2
# freqency scales
omega2=3.
ext_omega2=0.6

id=4

prepath='/Users/edoardo/Desktop/simulazione_prova/record/'
ID = str(P)+'-'+str(N) + '-' + str(dt) +'-'+ str(id) +'.txt'
path=prepath+ID

'''
Basically it's working kinda good, maybe it could be a good idea to create some automatic
test to see when it fails to chek the presence of minor problems.
'''


def initialize(positions = None):
    """
    return initialized data in shape (P,N,dim)
    still not tested for non random initialization
    more than 2 dimension not yet implemented
    """
    if positions == None:
        pos = np.random.rand(P,dim)*large_scale-large_scale/2
    else: 
        pos = position
    lattice = l*np.asarray([[np.cos(2*np.pi*i/N),np.sin(2*np.pi*i/N)] for i in range(N)])
    noise = np.random.rand(P,N,dim)*disp_scale-disp_scale/2
    v_n = np.random.rand(N,P,dim)*vel_beads_scale-vel_beads_scale/2
    q_n = (lattice + noise).swapaxes(0,1)+pos
    return q_n, v_n

def static_plot(q):
    x = q.swapaxes(0,2)[0].flatten()
    y = q.swapaxes(0,2)[1].flatten()
    fig, ax = plt.subplots()
    foo = large_scale
    ax.axis([-foo, foo, -foo, foo])
    ax.scatter(x,y,c='#2300A8')
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
    plt.show()

   
#static_plot(initialize())


def dV(q):
    """
    If you want the free case just set ext_omega2 == 0.
    """
    dv = np.zeros((N,P,dim))
    for i in range(N):
        dv[i] = omega2*(2*q[i]-q[(i-1)%N]-q[(i+1)%N])+ext_omega2*q[i]
        #dv[i]=ext_omega2*q[i]
    return dv
 
'''    
def verlet_step(q,v,dV0,fix=None,pos=None,fixc=None,posc=None):
    """
    at the moment it seem like it's analitically problematic to fix more than 
    one mass at the same time, I'm not sure how to solve this.
    Make an experiment: try to invert the matrix even if in theory you shouldn't
    Achtung: there are a few problem with n+-1/2 solved by adatting q, pay attention!
    """
    #dV0=dV(q)
    q1=q+v*dt-dV0/2*dt**2
    dV1=dV(q1)
    if fix:
        lamda_r = np.zeros((N,P,dim))
        lamda_v = np.zeros((N,P,dim))
        for i,j in fix:
            #pos=np.zeros((N,P,2))
            lamda_v[i,j]=2/dt*v[i,j]-dV0[i,j]
            lamda_r[i,j]=(q1[i,j]/dt+v[i,j]-dt*dV1[i,j])/(dt/2)-lamda_v[i,j]
            lamda_r[i,j] -= 2/dt**2*pos[i,j]
            #achtung to pos
    #elif fixc:
    #    for i in range(N):
    #        lamda_r[i][fixc] = 2/dt**2*(q.sum(0)[fixc]-posc+v.sum(0)[fixc]*dt-dt**2/2*dV0.sum(0)[fixc]        
    else:
        i,j = 0,0
        lamda_v=0.
        lamda_r=0.
    v12=v-dt/2*dV0-dt/2*lamda_r
    #v12[i,j]-=dt/2*lamda_r
    v1=v12-dt/2*dV1-dt/2*lamda_v
    #v1[i,j]-=dt/2*lamda_v
    return q1,v1,dV1
    
    
def verlet_algorithm(q_0,v_0,fix=None,pos=None):
    Q_n, V_n= [], []
    q_n, v_n = q_0, v_0
    dV0=dV(q_n)
    for _ in range(steps):
        Q_n.append(q_n)
        V_n.append(v_n)
        if fix:
            pos = np.zeros((N,P,dim))
            for i,j in fix:
                #i,j=fix[0],fix[1]
                pos[i,j] = q_0[i,j]
        q_n, v_n, dV0= verlet_step(q_n, v_n,dV0,fix, pos)
    return np.asarray(Q_n), np.asarray(V_n)
    

'''


def constrained_r(q,v,dV0,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None):
    '''
    linear constrain can be combined simply adding forces, i still don't know if that's
    true in general
    '''
    lamda_r = np.zeros((N,P,dim))
    # evaluate lamda_r(n)
    if fix:
        for i,j in fix:
            lamda_r[i,j]+=2/dt**2*(q[i,j]-pos[i,j])+2/dt*(v[i,j]-dt/2*dV0[i,j])
    if fixc:
        q_=q.swapaxes(0,1)
        v_=v.swapaxes(0,1)
        u = np.zeros((P,dim))
        vu = np.zeros((P,dim))
        for j in fixc:
            u[j]=q_[j].sum(0)
            vu[j]=v_[j].sum(0)
            for i in range(N):
                lamda_r[i,j]+=1/N*(2/dt**2*(u[j]-posc[j])+2/dt*vu[j]-dV0.sum(0)[j])
    if lenf:
        B = np.zeros((N,P,dim))
        C = np.zeros((N,P,dim))
        lamda_rx=np.zeros((N,P))
        for (i,a),(j,b),l in lenf:
            B[i,a]=q[i,a]-q[j,b]+dt*(v[i,a]-v[j,b])-dt**2/2*(dV0[i,a]-dV0[j,b])
            C[i,a]=-dt**2*(q[i,a]-q[j,b])
        BC=(B*C).sum(-1)
        B2=(B**2).sum(-1)
        C2=(C**2).sum(-1)
        lamda_rp=(-BC+np.sqrt(BC**2-C2*(B2-l**2)))/C2
        lamda_rm=(-BC-np.sqrt(BC**2-C2*(B2-l**2)))/C2
        for (i,a),(j,b),_ in lenf:
            if abs(lamda_rp[i,a])<abs(lamda_rm[i,a]):
                lamda_rx[i,a]=lamda_rp[i,a]
            else:
                lamda_rx[i,a]=lamda_rm[i,a]
            #I call the following lamda but actually they are G*lambda
            lamda_r[i,a]+=(q[i,a]-q[j,b])*lamda_rx[i,a]
            lamda_r[j,b]+=-(q[i,a]-q[j,b])*lamda_rx[i,a]
    if lenfc:
        #print('here')
        B = np.zeros((P,dim)) 
        C = np.zeros((P,dim))
        lamda_rx=np.zeros((P,))
        u = q.sum(0)/N
        vu = v.sum(0)/N 
        dV0u = dV0.sum(0)/N
        for a,b,l in lenfc:
            B[a]=u[a]-u[b]+dt*(vu[a]-vu[b])-dt**2/2*(dV0u[a]-dV0u[b])
            C[a]=-dt**2*(u[a]-u[b])/N
        BC=(B*C).sum(-1)
        B2=(B**2).sum(-1)
        C2=(C**2).sum(-1)
        lamda_rp=(-BC+np.sqrt(BC**2-C2*(B2-l**2)))/C2
        lamda_rm=(-BC-np.sqrt(BC**2-C2*(B2-l**2)))/C2
        for a,b,_ in lenfc:
            if abs(lamda_rp[a])<abs(lamda_rm[a]):
                lamda_rx[a]=lamda_rp[a]
            else:
                lamda_rx[a]=lamda_rm[a]
            #I call the following lamda but actually they are G*lambda
            for i in range(N):
                lamda_r[i,a]+=(u[a]-u[b])/N*lamda_rx[a]
                lamda_r[i,b]+=-(u[a]-u[b])/N*lamda_rx[a]
    return lamda_r

def constrained_v(q1,vp,dV1,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None):
    lamda_vp = np.zeros((N,P,dim))
    if fix:
        for i,j in fix:
            lamda_vp[i,j] += 2/dt*vp[i,j]-dV1[i,j]
    if fixc:
        vp_=vp.swapaxes(0,1)
        vpu=np.zeros((P,dim))
        for j in fixc:
            vpu[j]=vp_[j].sum(0)
            for i in range(N):
                lamda_vp[i,j]+=1/N*(2/dt*vpu[j]-dV1.sum(0)[j])
    if lenf:
        lamda_vx= np.zeros((N,P))
        for (i,a),(j,b),_ in lenf:
            q1vp=((q1[i,a]-q1[j,b])*(vp[i,a]-vp[j,b])).sum(-1)
            q1dV1=((q1[i,a]-q1[j,b])*(dV1[i,a]-dV1[j,b])).sum(-1)
            q2=((q1[i,a]-q1[j,b])**2).sum(-1)
            lamda_vx[i,a]=1/dt*(q1vp-dt/2*q1dV1)/(q2)
            lamda_vp[i,a]+=lamda_vx[i,a]*(q1[i,a]-q1[j,b])
            lamda_vp[j,b]+=-lamda_vx[i,a]*(q1[i,a]-q1[j,b])
    if lenfc:
        u1 = q1.sum(0)/N
        vpu = vp.sum(0)/N 
        dV1u = dV1.sum(0)/N
        lamda_vx= np.zeros((P,))
        #print(lenfc)
        for a,b,_ in lenfc:
            u1vpu=((u1[a]-u1[b])*(vpu[a]-vpu[b])).sum(-1)
            u1dV1u=((u1[a]-u1[b])*(dV1u[a]-dV1u[b])).sum(-1)
            u2=((u1[a]-u1[b])**2).sum(-1)
            lamda_vx[a]=1/dt*(u1vpu-dt/2*u1dV1u)/(u2)
            for i in range(N):
                lamda_vp[i,a]+=lamda_vx[a]*(u1[a]-u1[b])
                lamda_vp[i,b]+=-lamda_vx[a]*(u1[a]-u1[b])
    return lamda_vp
    
def verlet_step(q,v,dV0,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None):
    """
    vm == v_n-1/2
    vp == v_n+1/2
    lamda_v = lamda_v_n from the previous step
    lamda_vp = new lamda_v_(n+1)
    Achtung: there are a few problem with n+-1/2 solved by adatting q, pay attention!
    """
    # evaluate lamda_r(n)
    lamda_r = constrained_r(q,v,dV0,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1/2)
    vp = v-dV0*dt/2-dt/2*lamda_r
    # evaluate q(n+1)
    q1 = q+dt*vp
    dV1=dV(q1)
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_v(q1,vp,dV1,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1)
    v1=vp-dt/2*dV1-dt/2*lamda_vp
    return q1,v1,dV1
    
def verlet_algorithm(q_0,v_0,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None):
    Q_n, V_n= [], []
    q_n, v_n = q_0, v_0
    dV0=dV(q_n)
    for _ in range(steps):
        Q_n.append(q_n)
        V_n.append(v_n)
        if fix:
            pos = np.zeros((N,P,dim))
            for i,j in fix:
                pos[i,j] = q_0[i,j]   #probably not necessary
        if fixc:
            posc = np.zeros((P,dim))
            for j in fixc:
                posc[j]=q_0.swapaxes(0,1)[j].sum(0)
        q_n, v_n, dV0 = verlet_step(q_n, v_n,dV0,fix=fix,pos=pos,fixc=fixc,posc=posc,lenf=lenf,lenfc=lenfc)
    return np.asarray(Q_n), np.asarray(V_n)
    

def energy(Q_n,V_n):
    kin = (V_n**2).sum((1,2,3))/2
    pot1 = ext_omega2/2*(Q_n**2).sum((1,2,3))
    pot2 = 0.
    Q_n=Q_n.swapaxes(0,1)
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))
    pot2 *= omega2/2
    return kin+pot1+pot2
    
def angular_momentum(Q_n,V_n):
    return np.cross(Q_n,V_n).sum((1,2))
    
    
q_0,v_0 = initialize()

#fix=[(0,0),(3,0)]
#fixc=[0,2]

#lint=np.sqrt(((q_0[0,0]-q_0[1,0])**2).sum(-1))
#lenf=[((0,0),(1,0),lint)]
fix=None
#fixc=None
lenf=None
#lenfc=None

u_0=q_0.sum(0)
lintc=np.sqrt(((u_0[0]-u_0[1])**2).sum(-1))/N

lenfc=[(0,1,lintc)]
fixc =[0]

Q_n,V_n= verlet_algorithm(q_0,v_0,fix=fix,fixc=fixc,lenf=lenf,lenfc=lenfc)

G_n=(Q_n.swapaxes(0,1)).swapaxes(1,2)
#l_n = np.sqrt(((G_n[0,0]-G_n[1,0])**2).sum(-1))

E_n = energy(Q_n,V_n)
L_n = angular_momentum(Q_n,V_n)

out_file=open(path,'wb+')
pickle.dump([Q_n,V_n,E_n,L_n,dt],out_file)
out_file.close()

fig, (ax1,ax2) = plt.subplots(2,1)
x = np.linspace(0,steps*dt,steps-3)
E_n=E_n[3:]
L_n=L_n[3:]

#L_n=l_n[3:]-l_n.mean()
ax1.scatter(x,E_n,c='yellow',s=0.1)
ax2.scatter(x,L_n,c='red',s=0.1)
ax1.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax2.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.show()



'''
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
'''