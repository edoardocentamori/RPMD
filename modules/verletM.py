import numpy as np
from setting import *
from modules.rattle import constrained_r, constrained_v
from modules.normal import get_norm , prop_norm, get_stand, prop_normM
from modules.thermostat import thermal_step, debug_thermal_step, thermal_step2
from modules.verlet import *
from modules.rattleM import constrained_rM, constrained_vM
from modules.thermostatM import thermal_stepM, thermal_step2M

    
def verlet_stepM(q, v, dV0, beta, fix=None, pos=None, fixc=None, posc=None, lenf=None, lenfc=None):
    """
    vm == v_n-1/2
    vp == v_n+1/2
    lamda_v = lamda_v_n from the previous step
    lamda_vp = new lamda_v_(n+1)
    Achtung: there are a few problem with n+-1/2 solved by adatting q, pay attention!
    """
    # evaluate lamda_r(n)
    lamda_r = constrained_r(q, v, dV0, fix, pos, fixc, posc, lenf, lenfc)
    # evaluate v(n+1/2)
    vp = v-dV0/M*dt/2-dt/2*lamda_r
    # evaluate q(n+1)
    q1 = q+dt*vp
    dV1 = dV(q1)
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_v(q1, vp, dV1, fix, pos, fixc, posc, lenf, lenfc)
    # evaluate v(n+1)
    v1 = vp-dt/2*dV1/M-dt/2*lamda_vp
    return q1, v1, dV1

    
def verlet_step2M(q, v, dV0int, dV0ext, beta, fix=None, pos=None, fixc=None, posc=None, fixcm=None, poscm=None, lenf=None, lenfc=None, therm=None, H2=None):
    """
    New implementation of verlet_step using normal modes propagation.
    """
    # single thermalization in rattle
    if therm:
        vt = thermal_step2M(v, beta)
    else:
        vt = v
    # simulate thermal action as evolution through an effective potential    
    # dV_eff = -2*(vt-v)/dt
    # evaluate lamda_r(n) assuming the effective potential acted
    # modified with v-> vt in lambda_r
    #trial lamda -> lamda/M both for r and v
    #!!! achtung recently substituited lamda_r and v -> lamda_rM, it worked previously
    lamda_r = constrained_rM(q, vt, (dV0ext+dV0int)/M, fix, pos, fixc, posc, fixcm, poscm, lenf, lenfc)
    # evaluate v(n+1/2)    
    vp = vt-dV0ext/M*dt/2-dt/2*lamda_r/M
    # normal mode change variable
    u, vpu = get_norm(q, vp)
    # normal mode propagation
    u1, vpu1 = prop_normM(u, vpu)
    # back to standard coordinate
    q1, vp = get_stand(u1, vpu1)
    if H2:
        dV1ext, dV1int = dVH2(q1), dVint(q1)
    else:
        dV1ext, dV1int = dVext(q1), dVint(q1)
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_vM(q1, vp, (dV1ext+dV1int)/M, fix, pos, fixc, posc, fixcm, poscm, lenf, lenfc)
    # evaluate v(n+1)
    v1 = vp-dt/2*dV1ext/M-dt/2*lamda_vp/M
    return q1, v1, dV1int, dV1ext


def verlet_step1M(q, v, dV0ext, beta, fix=None, pos=None, fixc=None, posc=None, lenf=None, lenfc=None, therm=None):
    """
    new implementation of verlet_step using normal modes propagation.
    problem need to be fixed, try using python debugger
    """
    # first thermalization
    if therm:
        v = thermal_stepM(v, beta)
    # evaluate lamda_r(n)
    lamda_r = constrained_r(q, v, dV0ext, fix, pos, fixc, posc, lenf, lenfc)
    # evaluate v(n+1/2)    
    vp = v-dV0ext/M*dt/2-dt/2*lamda_r
    # normal mode change variable
    u, vpu = get_norm(q, vp)
    # normal mode propagation
    #need to implement M in here
    u1, vpu1 = prop_norm(u, vpu)
    # back to standard coordinate
    q1, vp = get_stand(u1, vpu1)
    dV1ext = dVext(q1)
    # evaluate lamda_v(n+1)
    lamda_vp = constrained_v(q1, vp, dV1ext, fix, pos, fixc, posc, lenf, lenfc)
    # evaluate v(n+1)
    v1 = vp-dt/2*dV1ext/M-dt/2*lamda_vp
    # second thermalization
    if therm:
        v1 = thermal_stepM(v1, beta)
    return q1, v1, dV1ext

    
def debug_verl_step(q, v, dV0ext, beta, fix=None, pos=None, fixc=None, posc=None, lenf=None, lenfc=None, therm=None):
    # In debug mode q and v are actually the normal modes
    #just to make make it work, useless
    dV1ext = dV0ext
    # first thermalization
    vp = debug_thermal_step(v)
    # normal mode propagation
    q1, v1 = prop_norm(q, vp)
    # second thermalization
    v1 = debug_thermal_step(v1)
    return q1, v1, dV1ext


def verlet_algorithmM(q_0, v_0, beta, fix=None, pos=None, fixc=None, posc=None, fixcm=None, poscm=None, lenf=None, lenfc=None, norm=0, therm=None, debug=False, H2=False):
    Q_n, V_n = [], []
    q_n, v_n = q_0, v_0
    dV0 = dV(q_n)
    if H2:
        dV0ext = dVH2(q_n)
    else: 
        dV0ext=dVext(q_n)
    dV0int = dVint(q_n)
    if debug:
        u_n, vu_n = get_norm(q_n, v_n)
    for _ in range(steps):
        Q_n.append(q_n)
        V_n.append(v_n)
        if fix:
            pos = np.zeros((N, P, dim))
            for i, j in fix:
                pos[i, j] = q_0[i, j]
        if fixc:
            posc = np.zeros((P, dim))
            for j in fixc:
                posc[j] = q_0.swapaxes(0,1)[j].sum(0)
        if fixcm:
            poscm = np.zeros(dim)
            mtot = 0.
            for a in fixcm:
                mtot += m[a]
                poscm += m[a]*q_0[:, a, :].sum(0)
            poscm /= mtot
        if debug: 
            u_n, vu_n, dV0 = debug_verl_step(u_n, vu_n, dV0, beta, fix=fix, pos=pos, fixc=fixc, posc=posc, lenf=lenf, lenfc=lenfc)
            q_n, v_n = get_stand(u_n, vu_n)
            continue
        if norm == 0:
            q_n, v_n, dV0 = verlet_stepM(q_n, v_n, dV0, beta, fix=fix, pos=pos, fixc=fixc, posc=posc, lenf=lenf, lenfc=lenfc)
        elif norm:
            # added dVint e verlet step 2
            q_n, v_n, dV0int,dV0ext = verlet_step2M(q_n, v_n, dV0int, dV0ext, beta, fix=fix, pos=pos, fixc=fixc, posc=posc, fixcm=fixcm, poscm=poscm, lenf=lenf, lenfc=lenfc, therm=therm, H2=H2)
    return np.asarray(Q_n), np.asarray(V_n)