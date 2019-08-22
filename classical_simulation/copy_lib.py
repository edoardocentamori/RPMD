import numpy as np
import matplotlib.pyplot as plt

# H2 bond strengh
fh2 = 1.
dt = 0.01
steps = 1000

# positions
q1 = np.ones(3)
q2 = np.ones(3)*5
Q = np.vstack((q1,q2))

# momenta

p1 = np.ones(3)*2
p2 = np.ones(3)
P = np.vstack((p1,p2))

m = np.asarray([1.,1.]) #hydrogen massess

M = np.broadcast_to(m.reshape(2,1),(2,3))

def cart_spher(xyz,pxyz):
    x,y,z = xyz #I'm using a non obvius casting when using array in here
    px,py,pz = pxyz 
    # position convertion
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arctan2(z, np.sqrt(x**2+y**2))
    phi = np.arctan2(y,x)
    # momenta convertion
    pr = np.cos(theta)*np.cos(phi)*px+np.cos(theta)*np.sin(phi)*py+np.sin(theta)*pz
    Ltheta = r*(-np.sin(theta)*np.cos(phi)*px-np.sin(theta)*np.sin(phi)*py+np.cos(theta)*pz)
    Lphi = r*np.cos(theta)*(-np.sin(phi)*px+np.cos(phi)*pz)
    return (r,theta,phi), (pr,Ltheta,Lphi)
    
def spher_cart(rthetaphi,Prthetaphi):
    r,theta,phi = rthetaphi
    pr,Lt,Lp = Prthetaphi
    Lt/=r
    Lp/=(r*np.cos(theta))
    # position convertion
    x = r*np.cos(theta)*np.cos(phi)
    y = r*np.cos(theta)*np.sin(phi)
    z = r*np.sin(theta)
    # momenta convertion
    px = np.cos(theta)*np.cos(phi)*pr - np.sin(theta)*np.cos(phi)*Lt-np.sin(phi)*Lp
    py = np.cos(theta)*np.sin(phi)*pr-np.sin(theta)*np.sin(phi)*Lt+np.cos(phi)*Lp
    pz = np.sin(theta)*pr+np.cos(theta)*Lt
    return (x,y,z),(px,py,pz)
    
def get_norm(Q,P):
    # handling position
    c = (m[0]*Q[0]+m[1]*Q[1])/(m[0]+m[1])
    d = Q[0]-Q[1]
    # handling momenta
    pc = P[0]+P[1]
    pd = m[1]/(m[0]+m[1])*P[0]-m[0]/(m[0]+m[1])*P[1]
    mid = cart_spher(d,pd)
    d, pd = d = np.asarray(mid[0]),np.asarray(mid[1])
    return np.vstack((c,d)),np.vstack((pc,pd))
    
def get_stand(N, PN):
    # unpacking values
    d,pd = np.asarray(spher_cart(N[1],PN[1]))
    c = N[0]
    pc = PN[0]
    # converting to 2body problem
    Q = np.vstack((c+m[1]/(m[0]+m[1])*d,c-m[0]/(m[0]+m[1])*d)) 
    P = np.vstack((m[0]/(m[0]+m[1])*pc+pd,(m[1]/(m[0]+m[1])*pc-pd)))
    return Q,P
    
def FH2(Q):
    return fh2*np.vstack(((Q[1]-Q[0]),(Q[0]-Q[1])))
 
def verlet_step(q,p,FH20):
    # evaluate lamda_r(n)
    #lamda_r = constrained_r(q,v,dV0,fix,pos,fixc,posc,lenf,lenfc)
    lamda_r = 0
    # evaluate v(n+1/2)
    p12 = p+FH20*dt/2-dt/2*lamda_r
    # evaluate q(n+1)
    q1 = q+dt*p12/M
    FH21=FH2(q1)
    # evaluate lamda_v(n+1)
    #lamda_vp = constrained_v(q1,vp,dV1,fix,pos,fixc,posc,lenf,lenfc)
    lamda_vp = 0
    # evaluate v(n+1)
    p1=p12+dt/2*FH21-dt/2*lamda_vp
    return q1,p1,FH21

def verlet_algorithm(q_0,p_0):
    Q_n, P_n= [], []
    q_n, p_n = q_0, p_0
    FH20 = FH2(q_0)
    for _ in range(steps):
        Q_n.append(q_n)
        P_n.append(p_n)
        q_n, p_n, FH20 = verlet_step(q_n, p_n, FH20)
    return np.asarray(Q_n), np.asarray(P_n)
    
def check_motion(Q_n):
    q1, q2 = Q_n[:,0,:], Q_n[:,1,:]
    l = np.sqrt(((q1-q2)**2).sum(-1))
    return l
    
def norm_analysis(Q_n,P_n):
    for q,p in Q_n,P_n
    
Q_n,P_n = verlet_algorithm(Q,P)

l = check_motion(Q_n)    
    
fig, ax = plt.subplots()

x = np.linspace(0,steps,steps)

ax.scatter(x,l,c='blue',s=0.1)
plt.show()
    