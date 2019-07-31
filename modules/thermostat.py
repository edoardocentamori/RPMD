from setting import *
from modules.normal import v_norm, v_stand, w
from numpy.random import normal
import numpy as np

gamma = 2*w
gamma[0] = 1/t0

c1=np.exp(-dt/2*gamma)
c2=np.sqrt(1-c1**2)

def thermal_step(v):
    vu = v_norm(v)
    vu = (c1*vu.swapaxes(0,-1)).swapaxes(0,-1)
    delta = (np.sqrt(1/beta)*c2*normal(size=(dim,P,N))).swapaxes(0,-1)
    vu+= delta
    return v_stand(vu) 