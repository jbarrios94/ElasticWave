# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:09:35 2021

@author: JuanBarrios
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import Graficar as gf

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

def Condicion_Inicial(NumVolx, NumVoly):
    q = np.zeros((3,NumVoly + 4,NumVolx + 4))
    return q

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

def Fimm2Dx(s, W1, W2):
    dim = np.shape(W1)
    F = np.zeros((dim[0], dim[1],dim[2] - 2))
    norm = np.sum(W1[:,:,1:-1]**2,0)
    norm += (norm == 0. ) * 1.
    theta = np.sum(W1[:,:,2:]*W1[:,:,1:-1],0)/norm
    W1L = Ultrabee(theta, CFL)*W1[:,:,1:-1]
    norm = np.sum(W2[:,:,1:-1]**2,0)
    norm += (norm == 0 ) * 1.
    theta = np.sum(W2[:,:,0:-2]*W2[:,:,1:-1],0)/norm
    W2L = Ultrabee(theta, CFL)*W2[:,:,1:-1]
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
    W1L = Ultrabee(theta, CFL)*W1[:,1:-1,:]
    norm = np.sum(W2[:,1:-1,:]**2,0)
    norm += (norm == 0 ) * 1.
    theta = np.sum(W2[:,0:-2,:]*W2[:,1:-1,:],0)/norm
    W2L = Ultrabee(theta, CFL)*W2[:,1:-1,:]
    F[:,...] = ((np.sign(-s[1:-2,:])*(1.-dt/dx*np.abs(-s[1:-2,:]))*W1L) +
        (np.sign(s[2:-1,:])*(1.-dt/dx*np.abs(s[2:-1,:]))*W2L))
    F = F * 0.5
    return F

def SuperBee(theta, cfl):
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

Ibc = 2 #tipo de condicion de frontera a la izquierda
Dbc = 2 #tipo de condicion de frontera a la derecha
Arbc = 2 #tipo de condicion de frontera a la derecha
Abbc = 3 #tipo de condicion de frontera a la derecha
"""enumeracion:
  1 condicion de frontera periodica
  2 condicion de frontera absorvente (zero-extrapolacion)
  3 condicion de rontera de pared solida (reflectivas)
"""

#-----------------------------------------------------------------------------
#Parametros Iniciales---------------------------------------------------------
#-----------------------------------------------------------------------------
rho_a = 1; rho_b = 1. #Densidad
k_a = 1; k_b = 1.  #Bulk modulus
Lim_Infx = 3000; Lim_Supx = 10000; NumVolx = 350 * 2**0 + 1
Lim_Infy = 0; Lim_Supy = 4200; NumVoly = 210 * 2**0 + 1 
t_inicial = 0; t_final = 1500
CFL = 0.45 #Condicion CFL para el calculo del paso de tiempo
Xc, Yc, Xn, Yn, dx, dy = Dominio(Lim_Infx, Lim_Supx, Lim_Infy, Lim_Supy, NumVolx, NumVoly) #Generar dominio computacional
c = np.ones((NumVolx + 4, NumVoly + 4)) * 1.5
# c[2:-2, 2:-2] = np.load("/home/juan/ElasticWave/Teste_seg_eage/seg_teste5_dx10.npy")
c = np.transpose(c)
c[0,:] = c[2,:]; c[-1,:] = c[-3,:]
c[1,:] = c[2,:]; c[-2,:] = c[-3,:]
c[:,0] = c[:,2]; c[:,-1] = c[:,-3]
c[:,1] = c[:,2]; c[:,-2] = c[:,-3]
b = 1 / (0.31 * (1e3*c)**0.25)
b[c < 1.51] = 1.0
Rho = 1/b
K = c**2 * Rho
#-----------------------------------------------------------------------------
VelMax = np.max(c) #Velocidad maxima para el calculo del CFL
dt = (dx * CFL)/VelMax #Tamaño del paso del tiempo
#dt = 2.316
q = Condicion_Inicial(NumVolx, NumVoly) #Calculo de la condicion inicial
q = BCx(q, Ibc, 0)
q = BCx(q, Dbc, 1)
q = BCy(q, Arbc, 1)
q = BCy(q, Abbc, 0)
qb = np.copy(q) #Copia de q, para evolucion temporal
Nt = int(round(t_final/dt))
receiver = np.zeros((Nt, NumVolx))
dim = np.shape(q)
dimx = (dim[0], dim[1],dim[2]-1)
dimy = (dim[0], dim[1]-1,dim[2])
isx = NumVolx//2
isy = NumVoly//3
#isy = np.where(np.linspace(0, 4200, NumVoly) == 40)[0][0]
time_values = np.arange(Nt)*dt
f0  = 0.005
t0 = 1./f0
a = 1.
r = (np.pi * f0 * (time_values - t0))
src = a * (1-2.*r**2)*np.exp(-r**2)

F = np.zeros((3,NumVoly+4,NumVolx+3))
G = np.zeros((3,NumVoly+3,NumVolx+4))

inicio = time.time()
for i in range(Nt):
    # t_inicial += dt
    W1, W3 = Riemann_Elasticity_2Dx(q, Rho, K, c)
    qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (W3[:,:,1:-2] + W1[:,:,2:-1])
    F[:,:,1:-1] = Fimm2Dx(c, W1, W3)
    qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (F[:,:,2:-1] - F[:,:,1:-2])
    qb = BCx(qb, Ibc, 0)
    qb = BCx(qb, Dbc, 1)
    W1, W3 = Riemann_Elasticity_2Dy(qb, Rho, K, c)
    qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (W3[:,1:-2,:] + W1[:,2:-1,:])
    G[:,1:-1,:] = Fimm2Dy(c, W1, W3)
    qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (G[:,2:-1,:] - G[:,1:-2,:])
    qb[0, isy + 2, isx + 2] = qb[0, isy + 2, isx + 2] + dt * src[i]/(K[isy + 2, isx + 2]*dx*dy)
    receiver[i, :] = (q[0]*K)[isy + 2,2:-2]
    qb = BCy(qb, Arbc, 1)
    qb = BCy(qb, Abbc, 0)
    q = np.copy(qb)
fin = time.time()
print (fin - inicio)

Stress = (q[0]*K)[2:-2,2:-2]
np.save("/home/juan/ElasticWave/Teste_seg_eage/p_dswpa_t1200_teste5_sb_dx10.npy", Stress)
np.save("/home/juan/ElasticWave/Teste_seg_eage/rec_dswpa_t1200_teste5_sb_dx10.npy", receiver)
gf.plot_image(Stress, extent = [0, Lim_Supx, Lim_Supy, 0])
is1 = np.where(np.linspace(3000, 10000, NumVolx) == 6600)[0][0]
gf.plot_seismic_traces([Stress[:, is1]], 0, t_final)
print(np.max(Stress[:, is1]))
plt.show()
