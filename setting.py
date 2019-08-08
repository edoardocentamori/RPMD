import numpy as np


# Energy conversion units

femtoSec2au = 41.341374578  # don't know
kcal2cm     = 349.75        # Kcal/mol -> cm^-1
au2cm       = 27.211*8065.5 # a.u -> cm^-1
ang2au      = 1.88973       # don't know
eV2cm       = 8065.6        # eV -> cm^-1


P = 1 # masses
N = 40 # beads per mass

steps = 100000
dim = 1

# time scales
dt=0.005
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

#omega2=3.
#ext_omega2=0.6

'''
hw = 100.
beta = 0.1
bhw = hw*beta

omega2 = N/bhw**2   #also equal to w^2/wext^2
ext_omega2 = 1.

'''

#w_e = 1. # remeber it's in femtosecond^-1
#ext_omega2 = w_e**2/femtoSec2au**2


ext_omega2=1.  #remove later
#ext_omega2=0.  #remove later

beta = 4 # It's in atomic units

#omega2 = N/beta**2*femtoSec2au**2

omega2 = N**2/beta**2 #remove later #added N -> N**2

'''
C'è la storia del beta_n, capire se può causare problemi
'''


# temperature scales


#T = 20. # I have the feeling 20. is the right order of magnitude
#beta = 1./T


id=16

# Energy conversion units

femtoSec2au = 41.341374578  # don't know
kcal2cm     = 349.75        # Kcal/mol -> cm^-1
au2cm       = 27.211*8065.5 # a.u -> cm^-1
ang2au      = 1.88973       # don't know
eV2cm       = 8065.6        # eV -> cm^-1


