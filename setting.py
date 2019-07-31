import numpy as np

P = 1 # masses
N = 6 # beads per mass

steps = 200000
dim = 2

# time scales
dt=0.0001
dt_res=1.
t0=0.1 #thermalizaton time for centroid, I have no idea of the range of values physically acceptable [but I guess 1-10 for now]

# lenght scales

large_scale=20. # scale in which masses are randomly sampled
disp_scale=0.8 # scale in which beads get moved from lattice randomly
l=3 # initial distance beads from centroid

# velocity scales

vel_beads_scale = 0.2

# freqency scales

omega2=3.
ext_omega2=0.6

# temperature scales


T = 20. # I have the feeling 20. is the right order of magnitude
beta = 1./T


id=6