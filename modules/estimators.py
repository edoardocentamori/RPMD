from setting import *
from modules.normal import w


def energy(Q_n, V_n):
    """
    Once I've seen it negative, there might be some problem with it
    """
    kin = (V_n**2).sum((1, 2, 3))/2
    pot1 = ext_omega2/2*(Q_n**2).sum((1, 2, 3))
    pot2 = 0.
    Q_n = Q_n.swapaxes(0, 1)
    for i in range(N):
        pot2 += ((Q_n[i]-Q_n[(i-1) % N])**2).sum((1, 2))
    pot2 *= omega2/2
    return kin+pot1+pot2


def angular_momentum(Q_n, V_n):
    return np.cross(Q_n, V_n).sum((1, 2))
    
 
def fictitious_energy(Q_n, V_n):
    kin = (V_n**2).sum((1, 2, 3))/2
    pot1 = ext_omega2/2*(Q_n**2).sum((1, 2, 3))
    pot2 = 0.
    Q_n=Q_n.swapaxes(0, 1)
    for i in range(N):
        pot2 += ((Q_n[i]-Q_n[(i-1) % N])**2).sum((1, 2))
    pot2 *= omega2/2
    return (kin+pot1+pot2)*hw

    
def primitive_energy(Q_n, V_n):
    """
    Pay attention for all estimators, not clear wheater it's necessary to 
    divide per N or per P*N in multidimensional system
    pot1 seem ok, pot2 seem too big
    It seem to mee that some particles are too disperded, maybe thermalization is wrong
    """
    pot1 = ext_omega2/2*(Q_n**2).sum((1, 2, 3))/N
    pot2 = 0.
    Q_n=Q_n.swapaxes(0, 1)
    for i in range(N):
        #pri=((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))[-1]
        #if pri >10:
        #    print(pri)
        pot2 += ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))
    pot2 *= omega2/(2*N) 
    #pot2 /= N #added as an experiment, need to check
    #print(pot1)
    #print('---')
    #print(pot2)
    #print('---')
    return N/(2*beta)-pot2+pot1


def real_primitive_energy(Q_n, V_n):
    pot1 = ext_omega2/2*(Q_n**2).sum((1, 2, 3))/N
    pot2 = 0.
    Q_n=Q_n.swapaxes(0, 1)
    for i in range(N):
        pot2 += ((Q_n[i]-Q_n[(i-1) % N])**2).sum((1, 2))
    pot2 *= omega2/(2*N) 
    kin = 1/2*(V_n**2).sum((1, 2, 3))/N
    # /N added to kin, since beads thermalize at beta/N not beta
    return kin-pot2+pot1


def real_primitive_energyH2(Q_n,V_n):
    pot1 = ext_omega2/2*((Q_n[:, :, 1, :]-Q_n[:, :, 0, :])**2).sum((2, 3))/N
    pot2 = 0.
    Q_n=Q_n.swapaxes(0, 1)
    for i in range(N):
        pot2 += ((Q_n[i]-Q_n[(i-1) % N])**2).sum((1, 2))
    pot2 *= omega2/(2*N) 
    kin= 1/2*(V_n**2).sum((1, 2, 3))/N
    # /N added to kin, since beads thermalize at beta/N not beta
    return kin-pot2+pot1


def virial_energy(Q_n, V_n):
    x_c = Q_n.sum(1)/N
    inter1 = (Q_n**2).sum((1, 2, 3))
    inter2 = ((Q_n.swapaxes(0, 1)-x_c).swapaxes(0, 1)*Q_n).sum((1, 2, 3))
    inter = (inter1+inter2)*ext_omega2
    return 1/beta*(inter/(2*N)+1/2.)
    

def mode_energy(U_n, VU_n):
    kin = (VU_n**2).sum((2, 3))/2
    pot = (U_n**2).sum((2, 3))/2*(w**2+ext_omega2)
    return kin+pot
    
    
def time_avarage(A):
    return A.cumsum()/(np.arange(len(A))+1.)
    
    
def kinetic(V_n):
    return (V_n**2).sum((1, 2, 3))/2/N
    
    
def pot1(Q_n):
    return ext_omega2/2*(Q_n**2).sum((1, 2, 3))/N
    
    
def pot2(Q_n):
    pot2 = 0.
    Q_n = Q_n.swapaxes(0, 1)
    for i in range(N):
        pot2 += ((Q_n[i]-Q_n[(i-1) % N])**2).sum((1, 2))
    pot2 *= omega2/2 
    return pot2/N