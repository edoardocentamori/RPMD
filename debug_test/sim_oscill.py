import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal, seed
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit

w = 3.162
steps = 100000
dt = 0.005

seed(2)

numerical_sin_w = 1/w*np.sin(w*dt)
   

def exact_step(q, v):
    v1 = np.cos(w*dt)*v-w*np.sin(w*dt)*q
    q1 = np.cos(w*dt)*q+numerical_sin_w*v
    return q1, v1
    
    
q_0, v_0 = 0., 0.
    
    
def motion(q_0, v_0, therm=False):
    q, v = q_0, v_0
    Q, V =[], []
    for i in range(steps):
        Q.append(q)
        V.append(v)
        if therm:
            v = thermal_step(v)
        q, v = exact_step(q, v) 
        if therm:
            v = thermal_step(v)
    return np.asarray(Q), np.asarray(V)
    
    
gamma = 2*w


c1 = np.exp(-dt/2*gamma)
c2 = np.sqrt(1-c1**2)

beta = 0.2


def thermal_step(v):
    delta = np.sqrt(1./beta)*c2*normal()
    return c1*v+delta
     

Q, V = motion(q_0, v_0, therm=True)

fig, ax = plt.subplots()

T = np.linspace(0, steps*dt, steps)

ax.scatter(T, Q, c='green', s=0.1)
ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.show()


def mode_energy(q, v):
    kin = v**2/2
    pot = q**2/2*w**2
    return kin+pot
    
    
E = mode_energy(Q, V)

bins = 20

fig, ax = plt.subplots()

ax.hist(E, bins, density=True, facecolor='b', alpha=0.75)

plt.show()

statistic, bin_edges, bin_number = binned_statistic(E, np.arange(bins), statistic='count', bins=bins)

# future warning: in future use of non-tuple sequence for multidimensional indexing will be deprecated

# statistic = binned statistic
# bin_edges = chosen edges
# bin_number = an array of size of the input array, each value is #the name of the bin corresponding to that element

y = statistic/len(E)
x = np.zeros(bins)

for i in range(len(bin_edges)-1):
    x[i] = (bin_edges[i]+bin_edges[i+1])/2


# funzione di fit

def F(x, C, b):
    return C*np.exp(-b*x)

    
init = (1., 1.)
popt, pcov = curve_fit(F, x, y, init)

print("C = ", popt[0], " +/- ", np.sqrt(pcov[0, 0]))
print("beta = ", popt[1], " +/- ", np.sqrt(pcov[1, 1]))

x1 = np.linspace(0, np.max(E), 10000)
y1 = F(x1, popt[0], popt[1])

fig1, ax1 = plt.subplots()

ax1.scatter(x, y, c='b')
ax1.scatter(x1, y1, c='g', s=0.01)
ax1.set_yscale('log')

plt.show()