import numpy as np

# Energy conversion units, use later to convert time in femptosecond and energy in something else

femtoSec2au = 41.341374578  # don't know
kcal2cm     = 349.75        # Kcal/mol -> cm^-1
au2cm       = 27.211*8065.5 # a.u -> cm^-1
ang2au      = 1.88973       # don't know
eV2cm       = 8065.6        # eV -> cm^-1


P = 2 # masses
N = 10 # beads per mass

steps = 100000
dim = 2

# time scales
dt=0.02
#dt_res=1.
t0=1. #thermalizaton time for centroid, I have no idea of the range of values physically acceptable 
#Adrien suggested 1/w_ext

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

ext_omega2=1.  #remove later

beta = 4 # It's in atomic units

omega2 = N**2/beta**2 #remove later

id=17



