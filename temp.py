# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Configuration step (Please run it before the simulation code!)

import numpy as np
import matplotlib
# Show Plot in The Notebook
#matplotlib.use("nbagg")
import matplotlib.pyplot as plt
import time

#matplotlib.rcParams['figure.facecolor'] = 'w'          # remove grey background

# Initialization of parameters

# Simple finite difference solver
# Elastic wave equation
# 1-D regular staggered grid

# Basic parameters
# nt = 1300                                              # number of time steps
# nx = 1000                                              # number of grid points in x
c0 = 2500                                              # velocity (m/sec) (shear wave)
eps = 0.5                                              # stability limit
# xmax = 1000000.                                        # maximum range (m)
rho0 = 2500.                                           # density (kg/m**3)
mu0 = rho0*c0**2.                                      # shear modulus (Pa)
nop = 4     

# Initialization of setup
# --------------------------------------------------------------------------
nx    = 401         # number of grid points 
# c0    = 2500        # acoustic velocity in m/s
rho    = 2500       # density in kg/m^3
Z0    = rho*c0      # impedance
# mu    = rho*c0**2   # shear modulus
# rho0  = rho         # density
# mu0   = mu          # shear modulus
isnap = 2                                              # snapshot frequency
isx = round(nx/2)                                      # source location
f0 = 1/15  
xmax  = 1000000.       # in m 
# eps   = 0.8         # CFL
tmax  = 50        # simulation time in s
# isnap = 10          # plotting rate
# sig   = 200         # argument in the inital condition
# x0    = 2500        # position of the initial condition     

# Finite Differences setup 
# --------------------------------------------------------------------------  
# dx    = xmax/(nx-1)                  # calculate space increment
# xfd   = np.arange(0, nx)*dx          # initialize space 
# mufd  = np.zeros(xfd.size) + mu0     # initialize shear modulus 
# rhofd = np.zeros(xfd.size) + rho0    # initialize density 

# Introduce inhomogeneity
#mufd[int((nx-1)/2) + 1:nx] = mufd[int((nx-1)/2) + 1:nx]*4

# initialize fields
# s  = np.zeros(xfd.size)
# v  = np.zeros(xfd.size)
# dv = np.zeros(xfd.size)
# ds = np.zeros(xfd.size)                                      # number of operator either 2 or 4
# 

dx = xmax / (nx - 1)                                     # calculate space increment (m)
x = (np.arange(nx)*dx)                                 # initialize space coordinates 
dt = eps * dx/c0

nt = round(tmax/dt)                                       # calculate time step from stability criterion(s)

# Source time function
# t = (np.arange(0,nt) * dt)                            # initialize time axis
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

# Extrapolation scheme and the plots

# Initialization of plot

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

# # Print setup parameters
# print("rho =", rho0, ", f_dom =", f0, ", stability limit =", eps, ", n_lambda", (lam/dx))

# # Initialize the plot
# title = "FD Elastic 1D staggered grid"
# fig = plt.figure(figsize=(12,8))
# ax1 = fig.add_subplot(2, 1, 1)
# ax2 = fig.add_subplot(2, 1, 2)
# line1 = ax1.plot(x, v, color = "red", lw = 1.5)
# line2 = ax2.plot(x, s, color = "blue", lw = 1.5)
# ax1.set_ylabel('velocity (m/s)')
# ax2.set_xlabel('x (m)')
# ax2.set_ylabel('stress (Pa)')
# plt.ion()
# plt.show()

# Begin extrapolation and update the plot
inicio = time.time()
for it in range (nt):
    
    # for i in range (2, nx-2):
    #     ds[i] = (0.0416666 * s[i-1] - 1.125 * s[i] + 1.125 * s[i+1] - 0.0416666 * s[i+2]) / (dx)
    
    # # Velocity extrapolation
    # v = v + dt * ds / rho
    
    # # Add source term at isx
    # v[isx] = v[isx] + dt * src[it] / (dt * rho[isx])
    
    # # Velocity derivative
    # for i in range (2, nx-2):

    #     dv[i] = (0.0416666 * v[i-2] - 1.125 * v[i-1] + 1.125 * v[i] - 0.0416666 * v[i+1]) / (dx)
        
    # Stress extrapolation
    # s = s + dt * mu * dv
        
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
    # ds[6:-6] = (2.18478116e-05 * s[1:-11] - 3.59005398e-04*s[2:-10] + 2.96728952e-03*s[3:-9] \
    #     -1.74476624e-02 * s[4:-8] + 9.69314575e-02 * s[5:-7] - 1.22133636e+00 * s[6:-6] \
    #         + 1.22133636e+00 * s[7:-5] -9.69314575e-02 * s[8:-4] + 1.74476624e-02 * s[9:-3] \
    #             -2.96728952e-03 * s[10:-2] + 3.59005398e-04 * s[11:-1] -2.18478116e-05 * s[12:])/ dx
    
    # v = v + dt * ds/rho
    
    # dv[6:-6] = (2.18478116e-05 * v[:-12] - 3.59005398e-04 * v[1:-11] + 2.96728952e-03*v[2:-10] \
    #     -1.74476624e-02 * v[3:-9] + 9.69314575e-02 * v[4:-8] - 1.22133636e+00 * v[5:-7] \
    #         + 1.22133636e+00 * v[6:-6] -9.69314575e-02 * v[7:-5] + 1.74476624e-02 * v[8:-4] \
    #             -2.96728952e-03 * v[9:-3] + 3.59005398e-04 * v[10:-2] -2.18478116e-05 * v[11:-1]) / dx
    
    # s = s + dt * mu * dv
    #------------------------------------------------------------------------
    
    
    
    #-------------------------------DF14-------------------------------------
    ds[7:-7] = (-4.23651475e-06 * s[1:-13] + 7.69225034e-05 * s[2:-12] - 6.89453549e-04 * s[3:-11] \
                +4.17893273e-03 * s[4:-10] - 2.04767704e-02 * s[5:-9] + 1.02383852e-01 * s[6:-8] \
                    -1.22860622e+00 * s[7:-7] + 1.22860622e+00 * s[8:-6] - 1.02383852e-01 * s[9:-5] \
                        +2.04767704e-02 * s[10:-4] -4.17893273e-03 * s[11:-3] + 6.89453549e-04 * s[12:-2] \
                            -7.69225034e-05 * s[13:-1] + 4.23651475e-06 * s[14:])/ dx
    
    v = v + dt * ds/rho
    
    dv[7:-7] = (-4.23651475e-06 * v[:-14] + 7.69225034e-05 * v[1:-13] - 6.89453549e-04 * v[2:-12] \
                +4.17893273e-03 * v[3:-11] - 2.04767704e-02 * v[4:-10] + 1.02383852e-01 * v[5:-9] \
                    -1.22860622e+00 * v[6:-8] + 1.22860622e+00 * v[7:-7] -1.02383852e-01 * v[8:-6] \
                        +2.04767704e-02 * v[9:-5] -4.17893273e-03 * v[10:-4] + 6.89453549e-04 * v[11:-3] \
                            -7.69225034e-05 * v[12:-2] + 4.23651475e-06 * v [13:-1]) / dx
    
    s = s + dt * mu * dv
    #------------------------------------------------------------------------
    
    
    # Stress extrapolation
    # s = s + dt * mu * dv
    
    s[isx] = s[isx] + dt * src[it]
    # 
    # Updating the plots
    # if not it % isnap: 
    #     for l in line1:
    #         l.remove()
    #         del l               
    #     for l in line2:
    #         l.remove()
    #         del l 
    #     line1 = ax1.plot(x, v, color = "red", lw = 1.5)
    #     line2 = ax2.plot(x, s, color = "blue", lw = 1.5)
    
    #     ax1.set_title(title + ", time step: %i" % (it))  
    #     plt.gcf().canvas.draw()
    
fin = time.time()

print(fin - inicio)

plt.plot(x, s)
plt.show()


