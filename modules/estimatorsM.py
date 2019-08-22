from setting import *
from modules.normal import w, W #not relly sure about W that has a strange definition

def energyM(Q_n,V_n):
    kin = (M*(V_n**2)).sum((1,2,3))/2
    pot1 = ext_omega2/2*(Q_n**2).sum((1,2,3))
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1)
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))
    pot2 *= omega2/2
    return kin+pot1+pot2
    
def angular_momentumM(Q_n,V_n):
    return (np.cross(Q_n,M*V_n)).sum((1,2))
    
 
def fictitious_energy(Q_n,V_n):
    kin = (M*(V_n**2)).sum((1,2,3))/2
    pot1 = ext_omega2/2*(Q_n**2).sum((1,2,3))
    pot2 = 0.
    Q_n=Q_n.swapaxes(0,1)
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))
    pot2 *= omega2/2
    return (kin+pot1+pot2)*hw
    
def primitive_energyM(Q_n,V_n):
    '''
    Pay attention for all estimators, not clear wheater it's necessary to 
    divide per N or per P*N in multidimensional system
    pot1 seem ok, pot2 seem too big
    It seem to mee that some particles are too disperded, maybe thermalization is wrong
    ACHTUNG : just added M on pot1
    '''
    pot1 = ext_omega2/2*(Q_n**2).sum((1,2,3))/N
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1) 
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2)) 
    pot2 *= omega2/(2*N) 
    return  dim*N/(2*beta)-pot2+pot1
    
def real_primitive_energyM(Q_n,V_n):
    '''
    I belive that pot2 should have M in it, but adding it trivially is not working !
    '''
    pot1 = ext_omega2/2*(Q_n**2).sum((1,2,3))/N
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1) 
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2)) 
    pot2 *= omega2/(2*N) 
    kin= (1/2*M*(V_n**2)).sum((1,2,3))/N
    # /N added to kin, since beads thermalize at beta/N not beta
    return kin-pot2+pot1
    
def real_primitive_energyH2M(Q_n,V_n):
    """
    ### added a /2 in pot1 is the same one coming from the reduced mass 
    """
    pot1 = ext_omega2/2/2*((Q_n[:,:,1,:]-Q_n[:,:,0,:])**2).sum((1,2))/N
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1) 
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2)) 
    pot2 *= omega2/(2*N) 
    kin= (1/2*M*(V_n**2)).sum((1,2,3))/N
    # /N added to kin, since beads thermalize at beta/N not beta
    return kin-pot2+pot1

def KH2M(Q_n,V_n):
    return (1/2*M*(V_n**2)).sum((1,2,3))/N

def pot1H2M(Q_n,V_n):
    return ext_omega2/2/2*((Q_n[:,:,1,:]-Q_n[:,:,0,:])**2).sum((1,2))/N

def pot2H2M(Q_n,V_n):
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1) 
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2)) 
    pot2 *= omega2/(2*N)
    return pot2


def virial_energyM(Q_n,V_n):
    x_c=Q_n.sum(1)/N
    inter1 = (Q_n**2).sum((1,2,3))
    inter2= ((Q_n.swapaxes(0,1)-x_c).swapaxes(0,1)*Q_n).sum((1,2,3))
    inter=(inter1+inter2)*ext_omega2
    # still not sure if 1/2 or N/2 are the correct things
    return  1/beta*(inter/(2*N)+1/2.)
    
def mode_energyM(U_n,VU_n):
    kin = (M*(VU_n**2)).sum((2,3))/2
    pot = ((M*(W**2)+ext_omega2)*U_n**2).sum((2,3))/2
    return kin+pot
    
def time_avarage(A):
    return A.cumsum()/(np.arange(len(A))+1.)
    
def kineticM(V_n):
    return (M*(V_n**2)).sum((1,2,3))/2/N
    
def pot1M(Q_n):
    return  ext_omega2/2*(Q_n**2).sum((1,2,3))/N
    
def pot2M(Q_n):
    pot2 = 0.
    Q_n=(np.sqrt(M)*Q_n).swapaxes(0,1)
    for i in range(N):
        pot2+= ((Q_n[i]-Q_n[(i-1)%N])**2).sum((1,2))
    pot2 *= omega2/2 
    return pot2/N