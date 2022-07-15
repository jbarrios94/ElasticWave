# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Configuration step (Please run it before the simulation code!)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time

c0 = 2.5                                               ## velocity (km/sec) (P-wave)
b = 1 / (0.31 * (1e3*c0)**0.25) 
rho0 = 1/b                                             
eps = 0.5                                              # stability limit
mu0 = rho0*c0**2.                                      # bulk modulus (K)   

# Initialization of setup
# --------------------------------------------------------------------------
nx    = 401                                            # number of grid points 
Z0    = rho0*c0                                        # impedance
isnap = 2                                              # snapshot frequency
isx = round(nx/2)                                      # source location
f0 = 1/15  
xmax  = 1000.                                          # in km 
tmax  = 100                                            # simulation time in s
dx = xmax / (nx - 1)                                   # calculate space increment (m)
x = (np.arange(nx)*dx)                                 # initialize space coordinates 
dt = eps * dx/c0                                       # calculate time step from stability criterion(s)
nt = round(tmax/dt)                                   

# Source time function
T0 = 1. / f0                                           # period
a = 4. / T0                                            # half-width (so called sigma)
t0 = T0 / dt
tmp = np.zeros(nt)
for it in range(nt):
    t = (it - t0) * dt
    tmp[it] = -2 * a * t * np.exp(-(a * t) ** 2)       # derivative of Gaussian (so called sigma)
src = np.zeros(nt)                                     # source
src[0:len(tmp)] = tmp
lam = c0 * T0                                          # wavelength

# Initialization of fields
v = np.zeros(nx)                                       # velocity
vnew = v
dv = v

s = np.zeros(nx)                                       # stress
snew = s
ds = s

mu = np.zeros(nx)                                      # shear modulus
rho = mu                                               
rho = rho + rho0
mu = mu + mu0

# Begin extrapolation and update the plot
inicio = time.time()
for it in range (nt):
    #-----------------------------DF4---------------------------------------
    # ds[2:-2] = (0.04166667 * s[1:-3] - 1.125 * s[2:-2] + 1.125* s[3:-1] - 0.04166667 * s[4:])/ dx
        
    # v = v + dt * ds/rho
    
    # dv[2:-2] = (0.04166667 * v[:-4] - 1.125 * v[1:-3] + 1.125* v[2:-2] - 0.04166667 * v[3:-1])/ dx
    
    # s = s + dt * mu * dv    
    #-----------------------------------------------------------------------
    
    
    
    #-----------------------------DF6---------------------------------------
    # ds[3:-3] = (-0.0046875 * s[1:-5] + 0.06510417 * s[2:-4] - 1.171875 * s[3:-3] \
    #     + 1.171875 * s[4:-2] - 0.06510417 * s[5:-1] + 0.0046875 * s [6:]) / dx
        
    # v = v + dt * ds/rho
    
    # dv[3:-3] = (-0.0046875 * v[:-6] + 0.06510417 * v[1:-5] - 1.171875 * v[2:-4] \
    #     + 1.171875 * v[3:-3] - 0.06510417 * v[4:-2] + 0.0046875 * v [5:-1]) / dx
    
    # s = s + dt * mu * dv    
    #-----------------------------------------------------------------------
    
    
    
    
    #------------------------------DF12 ------------------------------------
    ds[6:-6] = (2.18478116e-05 * s[1:-11] - 3.59005398e-04*s[2:-10] + 2.96728952e-03*s[3:-9] \
        -1.74476624e-02 * s[4:-8] + 9.69314575e-02 * s[5:-7] - 1.22133636e+00 * s[6:-6] \
            + 1.22133636e+00 * s[7:-5] -9.69314575e-02 * s[8:-4] + 1.74476624e-02 * s[9:-3] \
                -2.96728952e-03 * s[10:-2] + 3.59005398e-04 * s[11:-1] -2.18478116e-05 * s[12:])/ dx
    
    v = v + dt * ds/rho
    
    dv[6:-6] = (2.18478116e-05 * v[:-12] - 3.59005398e-04 * v[1:-11] + 2.96728952e-03*v[2:-10] \
        -1.74476624e-02 * v[3:-9] + 9.69314575e-02 * v[4:-8] - 1.22133636e+00 * v[5:-7] \
            + 1.22133636e+00 * v[6:-6] -9.69314575e-02 * v[7:-5] + 1.74476624e-02 * v[8:-4] \
                -2.96728952e-03 * v[9:-3] + 3.59005398e-04 * v[10:-2] -2.18478116e-05 * v[11:-1]) / dx
    
    s = s + dt * mu * dv
    #------------------------------------------------------------------------
    
    
    
    #-------------------------------DF14-------------------------------------
    # ds[7:-7] = (-4.23651475e-06 * s[1:-13] + 7.69225034e-05 * s[2:-12] - 6.89453549e-04 * s[3:-11] \
    #             +4.17893273e-03 * s[4:-10] - 2.04767704e-02 * s[5:-9] + 1.02383852e-01 * s[6:-8] \
    #                 -1.22860622e+00 * s[7:-7] + 1.22860622e+00 * s[8:-6] - 1.02383852e-01 * s[9:-5] \
    #                     +2.04767704e-02 * s[10:-4] -4.17893273e-03 * s[11:-3] + 6.89453549e-04 * s[12:-2] \
    #                         -7.69225034e-05 * s[13:-1] + 4.23651475e-06 * s[14:])/ dx
    
    # v = v + dt * ds/rho
    
    # dv[7:-7] = (-4.23651475e-06 * v[:-14] + 7.69225034e-05 * v[1:-13] - 6.89453549e-04 * v[2:-12] \
    #             +4.17893273e-03 * v[3:-11] - 2.04767704e-02 * v[4:-10] + 1.02383852e-01 * v[5:-9] \
    #                 -1.22860622e+00 * v[6:-8] + 1.22860622e+00 * v[7:-7] -1.02383852e-01 * v[8:-6] \
    #                     +2.04767704e-02 * v[9:-5] -4.17893273e-03 * v[10:-4] + 6.89453549e-04 * v[11:-3] \
    #                         -7.69225034e-05 * v[12:-2] + 4.23651475e-06 * v [13:-1]) / dx
    
    # s = s + dt * mu * dv
    #------------------------------------------------------------------------
    
    s[isx] = s[isx] + dt * src[it]

    
fin = time.time()

print(fin - inicio)
plt.plot(x, s)
plt.show()


