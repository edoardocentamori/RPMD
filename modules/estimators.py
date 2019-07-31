from setting import *


def energy(Q_n,V_n):
    '''
    Once I've seen it negative, there might be some problem with it
    '''
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