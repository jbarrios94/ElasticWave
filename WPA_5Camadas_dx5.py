# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:09:35 2021

@author: JuanBarrios
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

def DGauss(q, dt, t_inicial, Xc, Yc):
    a = 0.004
    f0 = 0.01
    src = -2.*a*(t_inicial - 1./f0)*np.exp(-a * (t_inicial - 1./f0)**2)
    c1 = Xc == 750.
    c2 = Yc == 750.
    q[0] = q[0] - (np.array([c1, c2]).all(0)*dt*src)
    return q

def Ricker(q, dt, t, isx, isz):
    t0 = 1./0.015
    a = 1.
    r = (np.pi * 0.015 * (t - t0))
    src = a * (1-2.*r**2)*np.exp(-r**2)
    q[0, isz, isx] = q[0, isx, isz] + dt*src
    return q

def Ricker(q, dt, t_inicial, Xc, Yc):
    t0 = 1./0.015
    # t0 = 0.
    a = 1.
    r = (np.pi * 0.015 * (t_inicial - t0))
    src = a * (1-2.*r**2)*np.exp(-r**2)
    c1 = Xc == Xc[3][np.shape(Xc)[1]//2]
    c2 = Yc == Yc[np.shape(Yc)[1]//2][3]
    q[0] = q[0] + (np.array([c1, c2]).all(0)*dt*src)
    return q

def Dominio(LimInfx, LimSupx, LimInfy, LimSupy, NumVolx, NumVoly):
    #Funcion para generar Dominio computacional
    #LimInf es el limite inferior
    #LimSup es el limite superior
    #NumVol es el numero de volumenes (celdas)
    dx = (LimSupx - LimInfx)/(NumVolx - 1)
    dy = (LimSupy - LimInfy)/(NumVoly - 1)
    xn = np.linspace(LimInfx, LimSupx + dx, NumVolx + 1)
    yn = np.linspace(LimInfy, LimSupy + dy, NumVoly + 1)
    xc = np.zeros(NumVolx + 4)
    yc = np.zeros(NumVoly + 4)
    xc[2:-2] = xn[:-1]
    yc[2:-2] = yn[:-1]
    Xc, Yc = np.meshgrid(xc,yc)
    Xn, Yn = np.meshgrid(xn, yn)
    return Xc, Yc, Xn, Yn, dx, dy

def Condicion_Inicial(Xc, Yc, Xn, Yn, NumVolx, NumVoly):
    q = np.zeros((3,NumVoly + 4,NumVolx + 4))
    return q

def velocidades(rho, k):
    c = np.sqrt(k/rho)
    return c

def Parametros_Dominio(Xc, Yc,rho_a, rho_b, k_a, k_b, NumVolx, NumVoly):
    #Teste1-------------------------------------------------------------------
    #Medio Heterogeneo--------------------------------------------------------
    # rho = rho_a * (Xc<=0) + rho_b * (Xc>0)
    # k = k_a * (Xc<=0) + k_b * (Xc>0)
    # Medio Homogeneo----------------------------------------------------------
    rho = rho_a * (Xc<=0) + rho_a * (Xc>0)
    k = k_a * (Xc<=0) + k_a * (Xc>0)
    #Teste2-------------------------------------------------------------------
    # rho = np.zeros(np.shape(Xc))
    # k = np.zeros(np.shape(Xc))
    # for j in range(NumVoly + 4):
    #     for i in range(NumVolx + 4):
    #         if i < np.max([(NumVolx + 4 )//2,j]) :   
    #             rho[j][i] = rho_a
    #             k[j][i] = k_a
    #         if i >= np.max([(NumVolx + 4)//2,j]) :   
    #             rho[j][i] = rho_b
    #             k[j][i] = k_b
    #         if j >= (NumVoly + 4)//2 and i==j :   
    #             rho[j][i] = (rho_a + rho_b)/2.
    #             k[j][i] = 1/(2*(k_a + k_b))
    #Medio Radial heterogeneo-------------------------------------------------
    # r = np.sqrt(Xc**2 + Yc**2)
    # a = (np.abs(r)> .5)
    # b = (np.abs(r)<= 1.)
    # c = (np.abs(r)<= .5)*(1.) + (np.array([a,b]).all(0))*(0.5) + (np.abs(r)> 1.)*(1.)
    # # c = (r<= 1.)*1. + (r>1.)*0.5
    # Z = c
    # # ## -----------------------------------------------------------------------------
    # rho = Z/c
    # k = Z * c
    return rho, k

def Riemann_Elasticity_2Dx(q, rho, k, c):
    r11 = 1/np.sqrt(k*rho)[:,:-1]
    r13 = -1/np.sqrt(k*rho)[:,1:]
    dfq1 = -((q[1]/rho)[:,1:] - (q[1]/rho)[:,:-1])
    dfq2 = -((q[0]*k)[:,1:] - (q[0]*k)[:,:-1])
    betha1 = (dfq1 - r13*dfq2)/(r11 - r13)
    betha3 = (-dfq1 + r11*dfq2)/(r11 - r13)
    W1 = np.zeros(dimx)
    W3 = np.zeros(dimx)
    W1[0] = betha1*r11
    W1[1] = betha1
    W3[0] = betha3*r13
    W3[1] = betha3
    return W1, W3


def Riemann_Elasticity_2Dy(q, rho, k, c):
    r11 = 1/np.sqrt(k*rho)[:-1,:]
    r13 = -1/np.sqrt(k*rho)[1:,:]
    dgq1 = -((q[2]/rho)[1:,:] - (q[2]/rho)[:-1,:])
    dgq3= -((q[0]*k)[1:,:] - (q[0]*k)[:-1,:])
    betha1 = (dgq1 - r13*dgq3)/(r11 - r13)
    betha3 = (-dgq1 + r11*dgq3)/(r11 - r13)
    W1 = np.zeros(dimy)
    W3 = np.zeros(dimy)
    W1[0] = betha1*r11
    W1[2] = betha1
    W3[0] = betha3*r13
    W3[2] = betha3
    return W1, W3

def Transvese_Riemann_Elasticity_2Dx(q, rho, k, Ap, Am, c):
    r11 = (1/np.sqrt(k*rho))[:-2,1:]
    r13 = (-1/np.sqrt(k*rho))[1:-1:,1:]
    gamma1 = (Ap[0,1:-1,:] - r13*Ap[2,1:-1,:])/(r11 - r13)
    dim = np.shape(r11)
    BmAp = np.zeros((3,dim[0],dim[1]))
    BmAp[0] = -c[:-2,1:]*gamma1*r11
    BmAp[2] = -c[:-2,1:]*gamma1

    r11 = (1/np.sqrt(k*rho))[1:-1,1:]
    r13 = (-1/np.sqrt(k*rho))[2:,1:]
    gamma3 = (-Ap[0,1:-1,:] + r11*Ap[2,1:-1,:])/(r11 - r13)
    dim = np.shape(r11)
    BpAp = np.zeros((3,dim[0],dim[1]))
    BpAp[0] = c[2:,1:]*gamma3*r13
    BpAp[2] = c[2:,1:]*gamma3
    
    r11 = (1/np.sqrt(k*rho))[:-2,:-1]
    r13 = (-1/np.sqrt(k*rho))[1:-1,:-1]
    gamma1 = (Am[0,1:-1,:] - r13*Am[2,1:-1,:])/(r11 - r13)
    dim = np.shape(r11)
    BmAm = np.zeros((3,dim[0],dim[1]))
    BmAm[0] = -c[:-2,:-1]*gamma1*r11
    BmAm[2] = -c[:-2,:-1]*gamma1
    
    r11 = (1/np.sqrt(k*rho))[1:-1,:-1]
    r13 = (-1/np.sqrt(k*rho))[2:,:-1]
    gamma3 = (-Am[0,1:-1,:] + r11*Am[2,1:-1,:])/(r11 - r13)
    dim = np.shape(r11)
    BpAm = np.zeros((3,dim[0],dim[1]))
    BpAm[0] = c[2:,:-1]*gamma3*r13
    BpAm[2] = c[2:,:-1]*gamma3
    return BpAp, BmAp, BpAm, BmAm


def Transvese_Riemann_Elasticity_2Dy(q, rho, k, Bp, Bm, c):
    r11 = (1/np.sqrt(k*rho))[1:,:-2]
    r13 = (-1/np.sqrt(k*rho))[1:,1:-1]
    gamma1 = (Bp[0,:,1:-1] - r13*Bp[1,:,1:-1])/(r11 - r13)
    dim = np.shape(r11)
    AmBp = np.zeros((3,dim[0],dim[1]))
    AmBp[0] = -c[1:,:-2]*gamma1*r11
    AmBp[1] = -c[1:,:-2]*gamma1
    
    r11 = (1/np.sqrt(k*rho))[1:,1:-1]
    r13 = (-1/np.sqrt(k*rho))[1:,2:]
    gamma3 = (-Bp[0,:,1:-1] + r11*Bp[1,:,1:-1])/(r11 - r13)
    dim = np.shape(r11)
    ApBp = np.zeros((3,dim[0],dim[1]))
    ApBp[0] = c[1:,2:]*gamma3*r13
    ApBp[1] = c[1:,2:]*gamma3
    
    r11 = (1/np.sqrt(k*rho))[:-1,:-2]
    r13 = (-1/np.sqrt(k*rho))[:-1,1:-1]
    gamma1 = (Bm[0,:,1:-1] - r13*Bm[1,:,1:-1])/(r11 - r13)
    dim = np.shape(r11)
    AmBm = np.zeros((3,dim[0],dim[1]))
    AmBm[0] = -c[:-1,:-2]*gamma1*r11
    AmBm[1] = -c[:-1,:-2]*gamma1
    
    r11 = (1/np.sqrt(k*rho))[:-1,1:-1]
    r13 = (-1/np.sqrt(k*rho))[:-1,2:]
    gamma3 = (-Bm[0,:,1:-1] + r11*Bm[1,:,1:-1])/(r11 - r13)
    dim = np.shape(r11)
    ApBm = np.zeros((3,dim[0],dim[1]))
    ApBm[0] = c[:-1,2:]*gamma3*r13
    ApBm[1] = c[:-1,2:]*gamma3
    return ApBp, AmBp, ApBm, AmBm

""""
all_white 
all_light_red 
all_light_blue 
all_light_green 
all_light_yellow 

red_white_blue 
blue_white_red 
red_yellow_blue 
blue_yellow_red 
yellow_red_blue 
white_red 
white_blue 

schlieren_grays 
schlieren_reds 
schlieren_blues
schlieren_greens
"""

def Fimm1D(s, W1, W2, L):
    norm = np.sum(W1[:,1:-1]**2,0)
    norm += (norm == 0. ) * 1.
    theta = np.sum(W1[:,2:]*W1[:,1:-1],0)/norm
    W1L = Limitador1D(theta, L)*W1[:,1:-1]
    norm = np.sum(W2[:,1:-1]**2,0)
    norm += (norm == 0 ) * 1.
    theta = np.sum(W2[:,0:-2]*W2[:,1:-1],0)/norm
    W2L = Limitador1D(theta, L)*W2[:,1:-1]
    F = ((np.sign(-s[1:-2])*(1.-dt/dx*np.abs(-s[1:-2]))*W1L) +
        (np.sign(s[2:-1])*(1.-dt/dx*np.abs(s[2:-1]))*W2L))
    F = F * 0.5
    return F

def Fimm2D(s, W1, W2, L, axe):
    dim = np.shape(W1)
    if axe == 1:
        F = np.zeros((dim[0], dim[1],dim[2] - 2))
        for j in range(np.shape(s)[0]):
            F[:,j,:] = Fimm1D(s[j,:], W1[:,j,:], W2[:,j,:], L)
    elif axe == 2:
        F = np.zeros((dim[0], dim[1] - 2,dim[2]))
        for j in range(np.shape(s)[1]):
            F[:,:,j] = Fimm1D(s[:,j], W1[:,:,j], W2[:,:,j], L)
    return F

def Fimm2Dx(s, W1, W2):
    dim = np.shape(W1)
    F = np.zeros((dim[0], dim[1],dim[2] - 2))
    norm = np.sum(W1[:,:,1:-1]**2,0)
    norm += (norm == 0. ) * 1.
    theta = np.sum(W1[:,:,2:]*W1[:,:,1:-1],0)/norm
    W1L = SuperBee(theta)*W1[:,:,1:-1]
    norm = np.sum(W2[:,:,1:-1]**2,0)
    norm += (norm == 0 ) * 1.
    theta = np.sum(W2[:,:,0:-2]*W2[:,:,1:-1],0)/norm
    W2L = SuperBee(theta)*W2[:,:,1:-1]
    F[:,...] = ((np.sign(-s[:,1:-2])*(1.-dt/dx*np.abs(-s[:,1:-2]))*W1L) +
        (np.sign(s[:,2:-1])*(1.-dt/dx*np.abs(s[:,2:-1]))*W2L))
    F = F * 0.5
    return F

def Fimm2Dy(s, W1, W2):
    dim = np.shape(W1)
    F = np.zeros((dim[0], dim[1] - 2,dim[2]))
    norm = np.sum(W1[:,1:-1,:]**2,0)
    norm += (norm == 0. ) * 1.
    theta = np.sum(W1[:,2:,:]*W1[:,1:-1,:],0)/norm
    W1L = SuperBee(theta)*W1[:,1:-1,:]
    norm = np.sum(W2[:,1:-1,:]**2,0)
    norm += (norm == 0 ) * 1.
    theta = np.sum(W2[:,0:-2,:]*W2[:,1:-1,:],0)/norm
    W2L = SuperBee(theta)*W2[:,1:-1,:]
    F[:,...] = ((np.sign(-s[1:-2,:])*(1.-dt/dx*np.abs(-s[1:-2,:]))*W1L) +
        (np.sign(s[2:-1,:])*(1.-dt/dx*np.abs(s[2:-1,:]))*W2L))
    F = F * 0.5
    return F

def Limitador1D(theta, o):
    shape = np.shape(theta)
    if o == 0:
        phy = 0.
        return phy
    if o == 1:
        phy = 1.
        return phy
    if o == 2:
        theta1 = np.array([np.ones(shape), theta])
        phy = MinMod(theta1)
        return phy
    if o == 3:
        a = np.zeros(shape)
        b = np.ones(shape)
        c = np.min([b , 2.*theta],0)
        d = np.min([2.*b, theta],0)
        phy = np.max([a,c,d],0)
        return phy
    if o == 4:
        a = np.zeros(shape)
        b = np.ones(shape)
        c = np.min([(b + theta)/2. , 2.*theta, 2*b],0)
        phy = np.max([a,c],0)
        return phy
    
def SuperBee(theta):
    shape = np.shape(theta)
    a = np.zeros(shape)
    b = np.ones(shape)
    c = np.min([b , 2.*theta],0)
    d = np.min([2.*b, theta],0)
    phy = np.max([a,c,d],0)
    return phy


def Ultrabee(theta, cfl):
    shape = np.shape(theta)
    ccfl = np.ones(shape) * cfl
    a = np.ones(shape) * 0.001
    b = np.ones(shape) * 0.999
    c = np.ones(shape)
    cfmod1 = np.max([a, ccfl], 0)
    cfmod2 = np.min([b, ccfl], 0)
    a = 2 * theta/cfmod1
    T1 = np.min([c, a], 0)
    b = 2/(1 - cfmod2)
    T2 = np.min([theta, b], 0)
    return np.max([np.zeros(shape), T1, T2], 0)
    

def cfl_superbee(r,cfl):
    r"""
    CFL-Superbee (Roe's Ultrabee) without theta parameter
    """

    a = np.empty((2,len(r)))
    b = np.zeros((3,len(r)))
    
    a[0,:] = 0.001
    a[1,:] = cfl
    cfmod1 = np.max(a[0],a[1])
    print(a.shape)
    print(cfmod1.shape)
    print(r.shape)
    a[0,:] = 0.999
    cfmod2 = np.min([a[0],a[1]], 1)
    
    a[0,:] = 1.0
    a[1,:] = 2.0 * r / cfmod1
    b[1,:] = np.min([a[0],a[1]], 0)
    a[0,:] = 2.0/(1-cfmod2)
    a[1,:] = r
    b[2,:] = np.min([a[0],a[1]], 0)
    
    return np.max(b,axis=0)

def MC(theta):
    shape = np.shape(theta)
    a = np.zeros(shape)
    b = np.ones(shape)
    c = np.min([(b + theta)/2. , 2.*theta, 2*b],0)
    phy = np.max([a,c],0)
    return phy

def MaxMod(a):
    k1 = a > 0.
    k2 = a < 0.
    return (k1.all(0))*np.max(a,0) + (k2.all(0))*np.min(a,0)

def MinMod(a):
    k1 = a > 0.
    k2 = a < 0.
    return (k1.all(0))*np.min(a,0) + (k2.all(0))*np.max(a,0)

def BCx(q, op, op1):
    if op1 == 0:
        if op == 1:
            q[:,:,0] = q[:,:,-4]
            q[:,:,1] = q[:,:,-3]
        if op == 2:
            q[:,:,0] = q[:,:,2]
            q[:,:,1] = q[:,:,2]
            return q
        if op == 3:
            q[:,:,1] = q[:,:,2]
            q[:,:,0] = q[:,:,3]
            q[1,:,0] = -q[1,:,3] 
            q[1,:,1] = -q[1,:,2]
            return q
    if op1 == 1:
        if op == 1:
            q[:,:,-1] = q[:,:,3]
            q[:,:,-2] = q[:,:,2]
        if op == 2:
            q[:,:,-1] = q[:,:,-3]
            q[:,:,-2] = q[:,:,-3]
            return q
        if op == 3:
            q[:,:,-2] = q[:,:,-3]
            q[:,:,-1] = q[:,:,-4]
            q[1,:,-1] = -q[1,:,-4] 
            q[1,:,-2] = -q[1,:,-3]
            return q

def BCy(q, op, op1):
    if op1 == 0:
        if op == 1:
            q[:,0,:] = q[:,-4,:]
            q[:,1,:] = q[:,-3,:]
        if op == 2:
            q[:,0,:] = q[:,2,:]
            q[:,1,:] = q[:,2,:]
            return q
        if op == 3:
            q[:,1,:] = q[:,2,:]
            q[:,0,:] = q[:,3,:]
            q[2,0,:] = -q[2,3,:] 
            q[2,1,:] = -q[2,2,:]
            return q
    if op1 == 1:
        if op == 1:
            q[:,-1,:] = q[:,3,:]
            q[:,-2,:] = q[:,2,:]
        if op == 2:
            q[:,-1,:] = q[:,-3,:]
            q[:,-2,:] = q[:,-3,:]
            return q
        if op == 3:
            q[:,-1,:] = q[:,-4,:]
            q[:,-2,:] = q[:,-3,:]
            q[2,-1,:] = -q[2,-4,:] 
            q[2,-2,:] = -q[2,-3,:]
            return q

def euler(des, vel, dt):
    des1 = des + dt * vel
    return des1



Ibc = 2 #tipo de condicion de frontera a la izquierda
Dbc = 2 #tipo de condicion de frontera a la derecha
Arbc = 2 #tipo de condicion de frontera a la derecha
Abbc = 2 #tipo de condicion de frontera a la derecha
"""enumeracion:
  1 condicion de frontera periodica
  2 condicion de frontera absorvente (zero-extrapolacion)
  3 condicion de rontera de pared solida (reflectivas)
"""
Lim = 3 #Limitador que sera usado para la reconstruccion
""" enumeracion
    0 primera orden (hace el termino de correccion F = 0)
    1 Metodo Lax-Wendroff
    2 Limitador de fluxo MinMod
    3 Limitador de fluxo SuperBee
    4 Limitador de fluxo MC-theta = 2
"""

#-----------------------------------------------------------------------------
#Parametros Iniciales---------------------------------------------------------
#-----------------------------------------------------------------------------
rho_a = 1; rho_b = 1. #Densidad
k_a = 1; k_b = 1.  #Bulk modulus
Lim_Infx = 0; Lim_Supx = 1000; NumVolx = 101
Lim_Infy = 0; Lim_Supy = 1000; NumVoly = 101
t_inicial = 0; t_final = 250
ricker = True
CFL = 0.45 #Condicion CFL para el calculo del paso de tiempo
Dimensional_Splitting = True; DS_ordem = 1
0
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
Xc, Yc, Xn, Yn, dx, dy = Dominio(Lim_Infx, Lim_Supx, Lim_Infy, Lim_Supy, NumVolx, NumVoly) #Generar dominio computacional
# Rho, K = Parametros_Dominio(Xc, Yc, rho_a, rho_b, k_a, k_b, NumVolx, NumVoly)
# c = velocidades(Rho, K) #Velocidades, son dadas al resolver el problema de riemann
#Ricker-----------------------------------------------------------------------
#Homogeneo--------------------------------------------------------------------
# c = np.ones(np.shape(Xc))*2.5
# Rho = np.ones(np.shape(Xc))
# K = c**2 * Rho
#-----------------------------------------------------------------------------
#Ricker-----------------------------------------------------------------------
# heterogeneo--------------------------------------------------------------------
# nlayers = 5
# cp_top = 1.5
# cp_bottom = 3.5
c = np.ones((NumVolx + 4, NumVoly + 4)) * 1.5
# c[:] = cp_top  # Top velocity (background)
# cp_i = np.linspace(cp_top, cp_bottom, nlayers)
# for i in range(1, nlayers):
#       c[i*int(NumVoly / nlayers) + 2:,:] = cp_i[i]  # Bottom velocity

# c[2:-2, 2:-2] = np.load('exemplo_Gaussianas.npy')
# c = np.transpose(c)
# c[0,:] = c[2,:]; c[-1,:] = c[-3,:]
# c[1,:] = c[2,:]; c[-2,:] = c[-3,:]
# c[:,0] = c[:,2]; c[:,-1] = c[:,-3]
# c[:,1] = c[:,2]; c[:,-2] = c[:,-3]


b = 1 / (0.31 * (1e3*c)**0.25)
b[c < 1.51] = 1.0
Rho = 1/b
K = c**2 * Rho
soluciones = []
#-----------------------------------------------------------------------------
VelMax = np.max(c) #Velocidad maxima para el calculo del CFL
dt = (dx * CFL)/VelMax #Tamaño del paso del tiempo
# dt = 1.6999999999999602
# dt = 0.00125
q = Condicion_Inicial(Xc, Yc, Xn, Yn, NumVolx, NumVoly) #Calculo de la condicion inicial
q = BCx(q, Ibc, 0)
q = BCx(q, Dbc, 1)
q = BCy(q, Arbc, 1)
q = BCy(q, Abbc, 0)
qb = np.copy(q) #Copia de q, para evolucion temporal
qini = np.copy(q)
desx = 0; desy = 0
Nt = int(round(t_final/dt))
# print(Nt)
# Nt = 0
receiver = np.zeros((Nt, NumVolx))
receiver1 = np.zeros((Nt, NumVolx))
receiver2 = np.zeros((Nt, NumVolx))


#Inicializacion de variables
dim = np.shape(q)
dimx = (dim[0], dim[1],dim[2]-1)
dimy = (dim[0], dim[1]-1,dim[2])
isx = round(NumVolx/2)
isy = round(NumVoly/2)

time_values = np.arange(Nt)*dt
f0  = 0.015
t0 = 1./f0
a = 1.
r = (np.pi * f0 * (time_values - t0))
src = a * (1-2.*r**2)*np.exp(-r**2)


F = np.zeros((3,NumVoly+4,NumVolx+3))
G = np.zeros((3,NumVoly+3,NumVolx+4))

if Dimensional_Splitting:
    if DS_ordem == 1:
        inicio = time.time()
        for i in range(Nt):
            # t_inicial += dt
            # print(t_inicial)
            W1, W3 = Riemann_Elasticity_2Dx(q, Rho, K, c)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (W3[:,:,1:-2] + W1[:,:,2:-1])
            F[:,:,1:-1] = Fimm2Dx(c, W1, W3)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (F[:,:,2:-1] - F[:,:,1:-2])
            # qb = BCx(qb, Ibc, 0)
            # qb = BCx(qb, Dbc, 1)
            W1, W3 = Riemann_Elasticity_2Dy(qb, Rho, K, c)
            qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (W3[:,1:-2,:] + W1[:,2:-1,:])
            G[:,1:-1,:] = Fimm2Dy(c, W1, W3)
            qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (G[:,2:-1,:] - G[:,1:-2,:])
            # qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2], dt, t_inicial, isx, isz)
            # qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2], dt, t_inicial, Xc[2:-2,2:-2], Yc[2:-2, 2:-2])
            qb[0, isy, isx] = qb[0, isy, isx] + dt * src[i]/K[isy, isx]
            # receiver[i, :] = q[0,3,2:-2]
            # receiver1[i, :] = q[1,3,2:-2]
            # receiver2[i, :] = q[2,3,2:-2]
            # if i == Nt//2:
            #Plot Stress-------------------------------------------------------
                # StressWPA = (q[0]*K)[2:-2,2:-2]
            # v_x = (q[1]/(Rho))[2:-2,2:-2]
            # v_z = (q[2]/(Rho))[2:-2,2:-2]
                # soluciones.append(StressWPA)
            # scale = .5*1e-3
            # plt.figure(figsize = (15, 15))
            # plt.imshow(StressWPA, cmap = 'RdGy', vmin = -10*scale, vmax = 10*scale)
            # plt.figure(figsize = (15, 15))
            # plt.imshow(v_x, cmap = 'RdGy', vmin = -scale, vmax = scale)
            # plt.figure(figsize = (15, 15))
            # plt.imshow(v_z, cmap = 'RdGy', vmin = -scale, vmax = scale)
            # qb = BCy(qb, Arbc, 0)
            # qb = BCy(qb, Abbc, 1)
            q = np.copy(qb)
        fin = time.time()
    else:
        for i in range(Nt):
            # print(t_inicial)
            t_inicial += dt
            F = np.zeros((3,NumVoly+4,NumVolx+3))
            G = np.zeros((3,NumVoly+3,NumVolx+4))
            W1, W3, Am, Ap = Riemann_Elasticity_2Dx(q, Rho, K, c)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/(2*dx) * (Ap[:,:,1:-2] + Am[:,:,2:-1])
            F[:,:,1:-1] = Fimm2D(c, W1, W3, Lim, 1)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/(2*dx) * (F[:,:,2:-1] - F[:,:,1:-2])
            qb = BCx(qb, Ibc, 0)
            qb = BCx(qb, Dbc, 1)
            W1, W3, Bm, Bp = Riemann_Elasticity_2Dy(qb, Rho, K, c)
            qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (Bp[:,1:-2,:] + Bm[:,2:-1,:])
            G[:,1:-1,:] = Fimm2D(c, W1, W3, Lim, 2)
            qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (G[:,2:-1,:] - G[:,1:-2,:])
            qb = BCy(qb, Arbc, 0)
            qb = BCy(qb, Abbc, 1)
            W1, W3, Am, Ap = Riemann_Elasticity_2Dx(qb, Rho, K, c)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/(2*dx) * (Ap[:,:,1:-2] + Am[:,:,2:-1])
            F[:,:,1:-1] = Fimm2D(c, W1, W3, Lim, 1)
            qb[:,:,2:-2] = qb[:,:,2:-2] - dt/(2*dx) * (F[:,:,2:-1] - F[:,:,1:-2])
            if ricker:
                qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2], dt, t_inicial, Xc[2:-2,2:-2], Yc[2:-2, 2:-2])
            receiver[i, :] = q[0,3,2:-2]
            receiver1[i, :] = q[1,3,2:-2]
            receiver2[i, :] = q[2,3,2:-2]
            #----------------------------------------------------------------------
            # if i == Nt//2:
            #Plot Stress-------------------------------------------------------
                # StressWPA = (q[0]*K)[2:-2,2:-2]
            # v_x = (q[1]/(Rho))[2:-2,2:-2]
            # v_z = (q[2]/(Rho))[2:-2,2:-2]
                # soluciones.append(StressWPA)
            # scale = .5*1e-3
            # plt.figure(figsize = (15, 15))
            # plt.imshow(StressWPA, cmap = 'RdGy', vmin = -10*scale, vmax = 10*scale)
            # plt.figure(figsize = (15, 15))
            # plt.imshow(v_x, cmap = 'RdGy', vmin = -scale, vmax = scale)
            # plt.figure(figsize = (15, 15))
            # plt.imshow(v_z, cmap = 'RdGy', vmin = -scale, vmax = scale)
            qb = BCx(qb, Ibc, 0)
            qb = BCx(qb, Dbc, 1)
            q = np.copy(qb)
else:
    for i in range(Nt):
        # if i%(Nt//10) == 0:
        #     np.save("FWPA_Stress_t"+str(round(t_inicial))+".npy",(q[0]*K)[2:-2,2:-2])
        #     # np.save("Vx1\FWPA_Vx_t"+str(round(t_inicial))+".npy",(q[1]/(Rho))[2:-2,2:-2])
        #     # np.save("Vz1\FWPA_Vz_t"+str(round(t_inicial))+".npy",(q[2]/(Rho))[2:-2,2:-2])
        #     # np.save("Receiver1\FWPA_Rec_t" + str(round(t_inicial)) + ".npy",receiver)
        # t_inicial += dt
        # print(t_inicial)
        F = np.zeros((3,NumVoly+4,NumVolx+3))
        G = np.zeros((3,NumVoly+3,NumVolx+4))
        W1x, W3x = Riemann_Elasticity_2Dx(q, Rho, K, c)
        F[:,:,1:-1] = F[:,:,1:-1] + Fimm2Dx(c, W1x, W3x)
        BpAp, BmAp, BpAm, BmAm = Transvese_Riemann_Elasticity_2Dx(q, Rho, K, W3x, W1x, c)
        G[:,:-1,1:] = G[:,:-1,1:] - dt/(2.*dx)*BmAp
        G[:,1:,1:] = G[:,1:,1:] - dt/(2.*dx)*BpAp
        G[:,:-1,:-1] = G[:,:-1,:-1] - dt/(2.*dx)*BmAm
        G[:,1:,:-1] = G[:,1:,:-1] - dt/(2.*dx)*BpAm
        
        W1y, W3y = Riemann_Elasticity_2Dy(q, Rho, K, c)
        G[:,1:-1,:] = G[:,1:-1,:] + Fimm2Dy(c, W1y, W3y)
        ApBp, AmBp, ApBm, AmBm = Transvese_Riemann_Elasticity_2Dy(q, Rho, K, W3y, W1y, c)
        F[:,1:,:-1] = F[:,1:,:-1] - dt/(2.*dy)*AmBp
        F[:,1:,1:] = F[:,1:,1:] - dt/(2.*dy)*ApBp
        F[:,:-1,:-1] = F[:,:-1,:-1] - dt/(2.*dy)*AmBm
        F[:,:-1,1:] = F[:,:-1,1:] - dt/(2.*dy)*ApBm
        
        
        qb[:,2:-2:,2:-2] = qb[:,2:-2,2:-2] - dt/dx * (W3x[:,2:-2,1:-2] + W1x[:,2:-2,2:-1])\
            - dt/dy * (W3y[:,1:-2,2:-2] + W1y[:,2:-1,2:-2])\
                - dt/dx * (F[:,2:-2,2:-1] - F[:,2:-2,1:-2])\
                    - dt/dy * (G[:,2:-1,2:-2] - G[:,1:-2,2:-2])\
                        
        # desx = euler(desx, (q[1]/(Rho))[2:-2,2:-2], dt)
        # desy = euler(desy, (q[2]/(Rho))[2:-2,2:-2], dt)
                        
        qb[0, isy, isx] = qb[0, isy, isx] + dt * src[i]                
        # qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2], dt, t_inicial, Xc[2:-2,2:-2], Yc[2:-2, 2:-2])
        # receiver[i, :] = q[0,3,2:-2]
        # qb = BCx(qb, Ibc, 0)
        # qb = BCx(qb, Dbc, 1)
        # qb = BCy(qb, Arbc, 1)
        # qb = BCy(qb, Abbc, 0)
        #----------------------------------------------------------------------
        q = np.copy(qb)
        

print (fin - inicio)
# np.save("Tiempo.txt" + str(dx),fin - inicio)
# np.save("FWPA_Stress_t"+str(round(t_inicial))+".npy",(q[0]*K)[2:-2,2:-2])
# np.save("Vx1\FWPA_Vx_t"+str(round(t_inicial))+".npy",(q[1]/(Rho))[2:-2,2:-2])
# np.save("Vz1\FWPA_Vz_t"+str(round(t_inicial))+".npy",(q[2]/(Rho))[2:-2,2:-2])
# np.save("Receiver1\FWPA_Rec_t" + str(round(t_inicial)) + ".npy",receiver)
#PColor Plot------------------------------------------------------------------
# fig, ax = plt.subplots()
# cnt = ax.pcolor(Xn, Yn, StressWPA, cmap = 'seismic')
# cnt = ax.imshow(StressWPA,vmin=-.5*1e-1, vmax=.5*1e-1, cmap = 'seismic')
# # plt.title('1')
# # ax.axvline(x=0, color = 'black')
# fig.colorbar(cnt)
# fig, ax = plt.subplots()
# # cnt = ax.pcolor(Xn, Yn, StressWPA, cmap = 'seismic')
# cnt = ax.imshow(v_x,vmin=-.5*1e-1, vmax=.5*1e-1, cmap = 'seismic')
# fig.colorbar(cnt)
# fig, ax = plt.subplots()
# # cnt = ax.pcolor(Xn, Yn, StressWPA, cmap = 'seismic')
# cnt = ax.imshow(v_z,vmin=-.5*1e-1, vmax=.5*1e-1, cmap = 'seismic')
# fig.colorbar(cnt)
# #------------------------------------------------------------------------------
# #Scatter plot-----------------------------------------------------------------
# R = np.reshape(np.sqrt(Xc[2:-2,2:-2]**2 + Yc[2:-2,2:-2]**2), NumVolx * NumVoly)
# S = np.reshape(StressWPA, NumVolx * NumVoly)
# # # print(np.abs(0.16056230868693167 - np.sum(np.abs(S))/(NumVolx * NumVoly)))
# fig1, ax1 = plt.subplots()
# cnt1 = ax1.plot(R, S, marker = '+', ms = 1.0, ls = "", label = "Dados 2D")
# cnt1 = ax1.pcolor(Xn, Yn, c[2:-2, 2:-2])
# plt.xlim(0., 3); plt.ylim(-1., 1.)
# plt.title('1')
#-----------------------------------------------------------------------------
#Schlieren Plot---------------------------------------------------------------
# (vx, vy) = np.gradient(StressWPA, dx, dy)
# vs = np.sqrt(vx**2 + vy**2)
# cnt = ax.contour(Xc[2:-2,2:-2], Yc[2:-2,2:-2], StressWPA, np.arange(-0.5,1.03,0.03), cmap = cm.schlieren_grays)
# # cnt = ax.pcolor(Xn, Yn, vs, cmap = cm.schlieren_grays)
# ax.plot([.0, .0],[-1, 0],'r', linewidth = 0.5)
# ax.plot([.0, 1],[0, 1],'r', linewidth = 0.5)
#plt.plot(xc[2:-2], Stress, 'r-', lw = 2.0, label = "Solução de referência 1D")
# print(format(np.sum(np.abs(S))/(NumVolx * NumVoly),'.3E'))
# Erro = np.abs(np.sum(np.abs(S))/(NumVolx* NumVoly))
# # print('N = ', NumVol)
# print('L1 = ',format(Erro, '.3E'))
# # Erro1 = np.abs(1.750000908615183 - np.max(Stress))
# Erro1 = np.abs(np.max(np.abs(S)))
# print('L-inf = ',format(Erro1, '.3E'))

#------------------------------------------------------------------------------
#Plot Densidade----------------------------------------------------------------
# plt.imshow(Rho[2:-2, 2:-2],cmap = 'jet',vmin=1., vmax = 3., extent = [0, 3000, 3000, 0])
#------------------------------------------------------------------------------
# Plot Stress-------------------------------------------------------------------
# scale = .5*1e-4
# plt.figure(figsize = (8, 8))
plt.imshow((q[0]*K)[2:-2,2:-2], cmap = 'seismic', extent = [0, Lim_Supx, Lim_Supy, 0])
plt.tight_layout()
plt.savefig("prueba3.png")

plt.figure()
plt.plot((q[0]*K)[NumVoly//2,2:-2], Xn[0,:-1])
plt.savefig("prueba4.png")
# plt.figure(figsize = (15, 15))
# plt.imshow(v_x, cmap = 'RdGy', vmin = -scale, vmax = scale)
# plt.figure(figsize = (15, 15))
# plt.imshow(v_z, cmap = 'RdGy', vmin = -scale, vmax = scale)
#Plot Receiver-----------------------------------------------------------------
# extent = [0., 70350., 1e-3*t_final, 0.]
# aspect = 9330./(1e-3*t_final)/.5
# plt.figure(figsize = (15, 15))
# plt.imshow(receiver, cmap = 'seismic', vmin = -.01, vmax = .01, interpolation = 'lanczos', extent=extent, aspect=aspect)
# print('Tiempo de ejecucion en segundos: ', fin-inicio)