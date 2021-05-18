import numpy as np

#This file is used to set the parameters of the system

# Energy conversion units, use later to convert time in femptosecond and energy in something else

femtoSec2au = 41.341374578  # femtosec to atomic units [a.u.]
kcal2cm     = 349.75        # Kcal/mol -> cm^-1
au2cm       = 27.211*8065.5 # a.u -> cm^-1
ang2au      = 1.88973       # don't know
eV2cm       = 8065.6        # eV -> cm^-1


m = [1.,1.] #set the masses of every beads
#m = [2000.,2000.] 
if len(m)==2:
    mu = m[0]*m[1]/(m[0]+m[1])
P = len(m) # masses
m = np.asarray(m)
N = 30 # beads per mass

steps = 1000
dim = 2

M = np.broadcast_to(m.reshape(1,P,1),(N,P,dim))


# time scales
dt=0.05 
ext_omega = 1.
#dt_res=1.
t0=1./ext_omega  # Thermalizaton time for centroid (Adrien suggested 1/w_ext)

# lenght scales

large_scale=20. # scale in which masses are randomly sampled
disp_scale=0.8 # scale in which beads get moved from lattice randomly
l=3 # initial distance beads from centroid

# velocity scales

vel_beads_scale = 0.2

# freqency scales
#w_e = 1. # remeber it's in femtosecond^-1
#ext_omega2 = w_e**2/femtoSec2au**2
#omega2 = N/beta**2*femtoSec2au**2

ext_omega2=ext_omega**2  

beta = 4 # It's in atomic units

omega2 = N**2/beta**2 #remove later

id=21 #identifier, unique used to differentiate between files with the same inputs.



