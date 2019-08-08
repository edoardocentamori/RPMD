import numpy as np
from setting import *
from modules.rattle import constrained_r, constrained_v
from modules.normal import get_norm , prop_norm, get_stand
from modules.thermostat import thermal_step, debug_thermal_step, thermal_step2



def zero_init():
    return np.zeros((N,P,dim)), np.zeros((N,P,dim))

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
    
def dV(q):
    """
    If you want the free case just set ext_omega2 == 0.
    """
    dv = np.zeros((N,P,dim))
    for i in range(N):
        dv[i] = omega2*(2*q[i]-q[(i-1)%N]-q[(i+1)%N])+ext_omega2*q[i]
    return dv
    
def dVext(q):
    """
    Just the external potential, probably easier to just call it directly
    """
    return ext_omega2*q
    
def verlet_step(q,v,dV0,beta,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None):
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
    
def verlet_step2(q,v,dV0ext,beta,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None, therm= None):
    """
    new implementation of verlet_step using normal modes propagation.
    problem need to be fixed, try using python debugger
    """
    # first thermalization
    if therm:
        v = thermal_step2(v,beta)
    # evaluate lamda_r(n)
    lamda_r = constrained_r(q,v,dV0ext,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1/2)    
    vp = v-dV0ext*dt/2-dt/2*lamda_r
    # normal mode change variable
    u, vpu = get_norm(q,vp)
    # normal mode propagation
    u1, vpu1 = prop_norm(u,vpu)
    # back to standard coordinate
    q1, vp = get_stand(u1,vpu1)
    dV1ext=dVext(q1)    
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_v(q1,vp,dV1ext,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1)
    v1=vp-dt/2*dV1ext-dt/2*lamda_vp
    return q1,v1,dV1ext
    
def verlet_step1(q,v,dV0ext,beta,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None, therm= None):
    """
    new implementation of verlet_step using normal modes propagation.
    problem need to be fixed, try using python debugger
    """
    # first thermalization
    if therm:
        v = thermal_step(v,beta)
    # evaluate lamda_r(n)
    lamda_r = constrained_r(q,v,dV0ext,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1/2)    
    vp = v-dV0ext*dt/2-dt/2*lamda_r
    # normal mode change variable
    u, vpu = get_norm(q,vp)
    # normal mode propagation
    u1, vpu1 = prop_norm(u,vpu)
    # back to standard coordinate
    q1, vp = get_stand(u1,vpu1)
    dV1ext=dVext(q1)    
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_v(q1,vp,dV1ext,fix,pos,fixc,posc,lenf,lenfc)
    # evaluate v(n+1)
    v1=vp-dt/2*dV1ext-dt/2*lamda_vp
    # second thermalization
    if therm:
        v1 = thermal_step(v1,beta)
    return q1,v1,dV1ext

    
def debug_verl_step(q,v,dV0ext,beta,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None, therm= None):
    # In debug mode q and v are actually the normal modes
    #just to make make it work, useless
    dV1ext=dV0ext
    # first thermalization
    vp = debug_thermal_step(v)
    # normal mode propagation
    q1, v1 = prop_norm(q,vp)
    # second thermalization
    v1 = debug_thermal_step(v1)
    return q1,v1,dV1ext
    
def verlet_algorithm(q_0,v_0,beta,fix=None,pos=None,fixc=None,posc=None,lenf=None,lenfc=None,norm=0,therm=None,debug=False):
    Q_n, V_n= [], []
    q_n, v_n = q_0, v_0
    dV0=dV(q_n) #problem in dimension different than 2 maybe fix later
    dV0ext=dVext(q_n)
    if debug:
        u_n, vu_n = get_norm(q_n,v_n)
    for _ in range(steps):
        Q_n.append(q_n)
        V_n.append(v_n)
        # remove this thing and add it in the constant
        if fix:
            pos = np.zeros((N,P,dim))
            for i,j in fix:
                pos[i,j] = q_0[i,j]   #probably not necessary
        if fixc:
            posc = np.zeros((P,dim))
            for j in fixc:
                posc[j]=q_0.swapaxes(0,1)[j].sum(0)
        if debug: 
            u_n, vu_n, dV0 = debug_verl_step(u_n, vu_n,dV0,beta,fix=fix,pos=pos,fixc=fixc,posc=posc,lenf=lenf,lenfc=lenfc)
            q_n,v_n = get_stand(u_n,vu_n)
            continue
        if norm==0:
            q_n, v_n, dV0 = verlet_step(q_n, v_n,dV0,beta,fix=fix,pos=pos,fixc=fixc,posc=posc,lenf=lenf,lenfc=lenfc)
        elif norm:
            q_n, v_n, dV0ext = verlet_step2(q_n, v_n,dV0ext,beta,fix=fix,pos=pos,fixc=fixc,posc=posc,lenf=lenf,lenfc=lenfc,therm=therm)
        #added verletstep2 to test it 
    return np.asarray(Q_n), np.asarray(V_n)