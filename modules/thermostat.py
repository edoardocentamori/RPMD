from setting import *
from modules.normal import v_norm, v_stand, w
from numpy.random import normal
import numpy as np
import matplotlib.pyplot as plt

gamma = 2*w
gamma[0] = 1/t0

c1 = np.exp(-dt/2*gamma)
c2 = np.sqrt(1-c1**2)


def thermal_step(v, beta):
    """
    It might have some problem for P>1 must fix
    """
    vu = v_norm(v)
    vu = (c1*vu.swapaxes(0, -1)).swapaxes(0, -1)
    ran = normal(size=(dim, P, N))
    delta = (np.sqrt(N/beta)*c2*ran).swapaxes(0, -1)
    vu += delta
    return v_stand(vu) 


def thermal_step2(v, beta):
    vu = v_norm(v)
    vu = (c1*vu.swapaxes(0, -1)).swapaxes(0, -1)
    vu += (np.sqrt(N/beta)*c2*normal(size=(dim, P, N))).swapaxes(0, -1)
    vu = (c1*vu.swapaxes(0, -1)).swapaxes(0, -1)
    vu += (np.sqrt(N/beta)*c2*normal(size=(dim, P, N))).swapaxes(0, -1)
    return v_stand(vu)

    
def debug_thermal_step(vu, beta):
    vu = (c1*vu.swapaxes(0, -1)).swapaxes(0, -1)
    ran = normal(size=(dim, P, N))
    delta = (np.sqrt(N/beta)*c2*ran).swapaxes(0, -1)
    vu += delta
    return vu 


fig, ax = plt.subplots()
x = np.linspace(0, N, N)

ax.scatter(x, c1, c='green', s=0.1, label='c1')
ax.scatter(x, c2, c='blue', s=0.1, label='c2')
ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax.set_xlabel('mode')
ax.set_ylabel('C1')
ax.legend()
ax.set_title('C1 and C2')
plt.show()