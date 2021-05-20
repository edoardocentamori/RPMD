import numpy as np
from setting import *
from modules.normal import S_,C_,alpha


def constrained_rM(q, v, dV0, fix=None, pos=None, fixc=None, posc=None, fixcm=None, poscm=None, lenf=None, lenfc=None):
    lamda_r = np.zeros((N, P, dim))
    if fix:  # fix singular beads
        for i, j in fix:
            lamda_r[i, j] += 2/dt**2*(q[i, j]-pos[i, j])+2/dt*(v[i, j]-dt/2*dV0[i, j])
    if fixc:  # fix centroids
        q_ = q.swapaxes(0, 1)
        v_ = v.swapaxes(0, 1)
        u = np.zeros((P, dim))
        vu = np.zeros((P, dim))
        for j in fixc:
            u[j] = q_[j].sum(0)
            vu[j] = v_[j].sum(0)
            # for i in range(N):
            #    lamda_r[i,j]+=1/N*(2/dt**2*(u[j]-posc[j])+2/dt*vu[j]-(dV0/M).sum(0)[j])
            # trial to test efficiency, not really more efficient but work
            lamda_r[:, j] += 1/N*(2/dt**2*(u[j]-posc[j])+2/dt*vu[j]-(dV0/M).sum(0)[j])
            # added dV0 -> dV0/M
    '''
    if fixcm:  # fixed center of mass
        mtot = 0.
        cm = np.zeros(dim)
        vcm = np.zeros(dim)
        dVm = np.zeros(dim)
        for a in fixcm:
            mtot += m[a]
            cm += m[a]*q[:,a,:].sum(0)
            vcm += m[a]*v[:,a,:].sum(0)
            dVm += dV0[:,a,:].sum(0)
        cm /= mtot
        vcm /= mtot
        for a in fixcm:
            lamda_r[:,a]+=(1/N*(2/dt**2*(cm-poscm)+2/dt*vcm-(dVm/mtot)))*mtot/m[a]
    '''
    if fixcm:  # new implementation, should work better in theory
        mtot = 0.
        q_ = np.einsum('ij,jlm->ilm', C_, q)
        v_ = np.einsum('ij,jlm->ilm', S_, v)
        dV_ = np.einsum('ij,jlm->ilm', S_, dV0)
        cm = np.zeros(dim)
        vcm = np.zeros(dim)
        dVm = np.zeros(dim)
        for a in fixcm:
            mtot += m[a]
            cm += m[a]*q_[:, a, :].sum(0)
            vcm += m[a]*v_[:, a, :].sum(0)
            dVm += dV_[:, a, :].sum(0)
        cm /= mtot
        vcm /= mtot
        dVm /= mtot
        for a in fixcm:
            lamda_r[:, a] += 2/(alpha*dt*N)*(cm-poscm+vcm-dVm*dt/2)*mtot/m[a]
    
    if lenf:  # fix lenght between beads
        B = np.zeros((N, P, dim))
        C = np.zeros((N, P, dim))
        lamda_rx = np.zeros((N, P))
        for (i, a), (j, b), l in lenf:
            B[i, a] = q[i, a]-q[j, b]+dt*(v[i, a]-v[j, b])-dt**2/2*(dV0[i, a]-dV0[j, b])
            C[i, a] = -dt**2*(q[i, a]-q[j, b])
        BC = (B*C).sum(-1)
        B2 = (B**2).sum(-1)
        C2 = (C**2).sum(-1)
        lamda_rp = (-BC+np.sqrt(BC**2-C2*(B2-l**2)))/C2
        lamda_rm = (-BC-np.sqrt(BC**2-C2*(B2-l**2)))/C2
        for (i, a), (j, b), _ in lenf:
            if abs(lamda_rp[i, a]) < abs(lamda_rm[i, a]):
                lamda_rx[i,a] = lamda_rp[i, a]
            else:
                lamda_rx[i, a] = lamda_rm[i, a]
            # I call the following lamda but actually they are G*lambda
            lamda_r[i, a] += (q[i, a]-q[j, b])*lamda_rx[i, a]
            lamda_r[j,b]+=-(q[i,a]-q[j,b])*lamda_rx[i, a]
    if lenfc:  # fix lenght between centroids
        # I'll try just dV -> dV/M without caring to details
        B = np.zeros((P, dim))
        C = np.zeros((P, dim))
        lamda_rx = np.zeros((P,))
        u = q.sum(0)/N
        vu = v.sum(0)/N 
        dV0u = (dV0/M).sum(0)/N # modified for M
        for a, b, l in lenfc:
            B[a] = u[a]-u[b]+dt*(vu[a]-vu[b])-dt**2/2*(dV0u[a]-dV0u[b])
            C[a] = -dt**2*(u[a]-u[b])/N
        BC = (B*C).sum(-1)
        B2 = (B**2).sum(-1)
        C2 = (C**2).sum(-1)
        lamda_rp = (-BC+np.sqrt(BC**2-C2*(B2-l**2)))/C2
        lamda_rm = (-BC-np.sqrt(BC**2-C2*(B2-l**2)))/C2
        for a, b, _ in lenfc:
            if abs(lamda_rp[a]) < abs(lamda_rm[a]):
                lamda_rx[a] = lamda_rp[a]
            else:
                lamda_rx[a] = lamda_rm[a]
            # I call the following lamda but actually they are G*lambda
            for i in range(N):
                lamda_r[i, a] += (u[a]-u[b])/N*lamda_rx[a]
                lamda_r[i, b] += -(u[a]-u[b])/N*lamda_rx[a]
    return lamda_r


def constrained_vM(q1, vp, dV1, fix=None, pos=None, fixc=None, posc=None, fixcm=None, poscm=None, lenf=None, lenfc=None):
    lamda_vp = np.zeros((N, P, dim))
    if fix:  # fix singular beads
        for i, j in fix:
            lamda_vp[i, j] += 2/dt*vp[i, j]-dV1[i, j]
    if fixc:  # fix centroids
        vp_ = vp.swapaxes(0, 1)
        vpu = np.zeros((P, dim))
        for j in fixc:
            vpu[j] = vp_[j].sum(0)
            for i in range(N):
                lamda_vp[i, j] += 1/N*(2/dt*vpu[j]-(dV1/M).sum(0)[j])
    if fixcm:
        mtot = 0.
        vcm = np.zeros(dim)
        dVm = np.zeros(dim)
        for a in fixcm:
            mtot += m[a]
            vcm += m[a]*vp[:, a, :].sum(0)
            dVm += dV1[:, a, :].sum(0)
        vcm /= mtot
        for a in fixcm:
            lamda_vp[:, a] += 1/N*(2/dt*vcm-(dVm/mtot))*mtot/m[a]
    if lenf:  # fix lenght between beads
        lamda_vx = np.zeros((N, P))
        for (i, a), (j, b), _ in lenf:
            q1vp = ((q1[i, a]-q1[j, b])*(vp[i, a]-vp[j, b])).sum(-1)
            q1dV1 = ((q1[i, a]-q1[j, b])*(dV1[i, a]-dV1[j, b])).sum(-1)
            q2=((q1[i, a]-q1[j, b])**2).sum(-1)
            lamda_vx[i, a] = 1/dt*(q1vp-dt/2*q1dV1)/q2
            lamda_vp[i, a] += lamda_vx[i, a]*(q1[i, a]-q1[j, b])
            lamda_vp[j, b] += -lamda_vx[i, a]*(q1[i, a]-q1[j, b])
    if lenfc:  # fix lenght between centroids
        u1 = q1.sum(0)/N
        vpu = vp.sum(0)/N 
        dV1u = (dV1/M).sum(0)/N  # modified for M
        lamda_vx = np.zeros((P,))
        for a, b, _ in lenfc:
            u1vpu = ((u1[a]-u1[b])*(vpu[a]-vpu[b])).sum(-1)
            u1dV1u = ((u1[a]-u1[b])*(dV1u[a]-dV1u[b])).sum(-1)
            u2 = ((u1[a]-u1[b])**2).sum(-1)
            lamda_vx[a] = 1/dt*(u1vpu-dt/2*u1dV1u)/u2
            for i in range(N):
                lamda_vp[i, a] += lamda_vx[a]*(u1[a]-u1[b])
                lamda_vp[i, b] += -lamda_vx[a]*(u1[a]-u1[b])
    return lamda_vp