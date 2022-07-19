#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:13:23 2022

@author: cbarrios
"""
r"""
Two-dimensional acoustics
=========================
Solve the (linear) acoustics equations:
.. math:: 
    p_t + K (u_x + v_y) & = 0 \\ 
    u_t + p_x / \rho & = 0 \\
    v_t + p_y / \rho & = 0.
Here p is the pressure, (u,v) is the velocity, K is the bulk modulus,
and :math:`\rho` is the density.
"""
from clawpack import riemann
import numpy as np
import matplotlib.pyplot as plt

from clawpack import pyclaw

def qinit(state,width=0.2):
    state.q[0,:,:] = 0.
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    
def ricker(solver,state,dt):
    """
    Geometric source terms for Euler equations with cylindrical symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    q = state.q
    f0  = 0.015
    t0 = 1./f0
    a = 1.
    r = (np.pi * f0 * (claw.solution.t - t0))
    src = a * (1-2.*r**2)*np.exp(-r**2)
    # print(src)
    q[0, 100, 100] = q[0, 100, 100] + dt * src
    
def dricker(solver,state,dt):
    """
    Geometric source terms for Euler equations with cylindrical symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    q = state.q
    f0  = 0.015
    t0 = 1./f0
    a = 1.
    r = (np.pi * f0 * (claw.solution.t - t0))
    src = a * (1-2.*r**2)*np.exp(-r**2)
    # print(src)
    dq = np.zeros(q.shape)
    dq[0, 100, 100] =  dt * src
    return dq
    

solver_type = 'sharpclaw'
time_integrator='SSP104'
if solver_type == 'classic':
    solver = pyclaw.ClawSolver2D(riemann.acoustics_2D)
    solver.dimensional_split=True
    solver.dt_variable = True
    solver.cfl_max = 1.
    solver.step_source = ricker
    # solver.cfl_max = 0.25
    solver.cfl_desired = 0.45
    # solver.limiters = 5
    solver.limiters = pyclaw.limiters.tvd.MC
elif solver_type=='sharpclaw':
    solver=pyclaw.SharpClawSolver2D(riemann.acoustics_2D)
    solver.time_integrator=time_integrator
    if solver.time_integrator=='SSP104':
        solver.dq_src = dricker
        solver.cfl_max = 0.5
        solver.cfl_desired = 0.45
    elif solver.time_integrator=='SSPLMMk2':
        solver.lmm_steps = 3
        solver.lim_type = 2
        solver.cfl_max = 0.25
        solver.cfl_desired = 0.24
    else:
        raise Exception('CFL desired and CFL max have not been provided for the particular time integrator.')


solver.bc_lower[0]=pyclaw.BC.extrap
solver.bc_upper[0]=pyclaw.BC.extrap
solver.bc_lower[1]=pyclaw.BC.extrap
solver.bc_upper[1]=pyclaw.BC.extrap

mx=200; my=200
x = pyclaw.Dimension(0.,1000.,mx,name='x')
y = pyclaw.Dimension(0.,1000.,my,name='y')
domain = pyclaw.Domain([x,y])

num_eqn = 3
state = pyclaw.State(domain,num_eqn)

# rho  = 1.0  # Material density
# bulk = 4.0  # Material bulk modulus
# cc = np.sqrt(bulk/rho)  # sound speed
cc = 2.5
rho = 1.
bulk = cc**2 * rho
zz = rho*cc             # impedance
state.problem_data['rho']= rho
state.problem_data['bulk']=bulk
state.problem_data['zz']= zz
state.problem_data['cc']=cc

solver.dt_initial=np.min(domain.grid.delta)/state.problem_data['cc']*solver.cfl_max

qinit(state)


claw = pyclaw.Controller()
claw.solution = pyclaw.Solution(state,domain)
claw.solver = solver
claw.output_format = None
claw.tfinal = 150
claw.run()

plt.imshow(state.q[0], cmap = 'seismic')



