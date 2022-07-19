# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:09:35 2021

@author: JuanBarrios
"""

#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
def ejecutar_dswpa(ref_factor):
    import time
    import numpy             as np
    import matplotlib.pyplot as plt
    #==============================================================================
    def Ricker(q, f0, dt, t_inicial, Xc, Yc, dx, dy, factor, yrec):
        t0 = 1./f0
        # t0 = 0.
        a = 1.
        r = (np.pi * f0 * (t_inicial - t0))
        src = a * (1-2.*r**2)*np.exp(-r**2)
        q[0] = q[0] + np.array([Yc==yrec,Xc== (Xc[0, Xc.shape[1]//2])]).all(0) * (4**factor*dt*src/K[2:-2,2:-2])
        return q
    
    def test3(ref_factor):
        
        rf = ref_factor
        
        shape      = (int(300* 2**rf + 1), int(100 * 2**rf + 1))
        spacing    = (10,10)
        origin     = (0,0)
        
        v          = np.empty(shape, dtype=np.float32)
        
        vel_media = np.array([1.5, 2.14648, 2.7168, 3.2056, 2.0654, 3.6873121, 4.6545323, 6.546846])
        
        layer_type = np.array([0, 1, 0, 0, 2, 2, 1, 1]);
        
        n_layers = 5
        
        prof = np.array([3, 3, 10, 5, -18, 0, 10, -10]) * 2**rf
        
        for k in range(n_layers):
            
            prof[k] = prof[k] + np.rint(shape[1]/n_layers)*k;
        
        k_vec = np.arange(0,n_layers)
        i_vec = np.arange(0,shape[0])
        j_vec = np.arange(0,shape[1])
        
        eq_layer = np.zeros((n_layers, shape[0]))
        
        incl = np.array([0.,0.05,-0.1,0.0,-0.0,-0.0, -0.1, 0.1])
        ampl = np.array([0,20,2,10,30,18,10,10]) * 2**rf
        wlen = np.array([0,0.2,1,0.8,0.06,0.3,0.5, 0.3])/2**rf
        cent = np.array([0, shape[0]*(100/300), shape[0]/4, shape[0]*(250/300),shape[0]/2, shape[0]/2, shape[0]*(300/300) , 0])
    
        for k in range(0,n_layers):
            
            if(layer_type[k]==0):
            
                eq_layer[k,:] = incl[k] * i_vec + prof[k];
            
            elif(layer_type[k]==1):
            
                eq_layer[k,:] = -ampl[k] * np.exp(-0.005/2**rf * wlen[k] * (i_vec-cent[k]) * (i_vec-cent[k])) + prof[k] + incl[k] * i_vec;
            else:
                eq_layer[k,:] = -ampl[k] * -np.exp(-0.005/2**rf * wlen[k] * (i_vec-cent[k]) * (i_vec-cent[k])) + prof[k] + incl[k] * i_vec;
        
        i_vec = np.array([i_vec,]*shape[1]).transpose()
        j_vec = np.array([j_vec,]*shape[0])
        
        v[:, :] = vel_media[0]
        
        for k in range(0, n_layers):
        
            indices = j_vec >= eq_layer[k, i_vec]
            
            v[np.where(indices)] = vel_media[k]
            
        return v
    
    
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
    
    def Riemann_Elasticity_2Dx(q, rho, k, c):
        r11 = 1/np.sqrt(k*rho)[:,:-1]
        r13 = -1/np.sqrt(k*rho)[:,1:]
        dfq1 = -((q[1]/rho)[:,1:] - (q[1]/rho)[:,:-1])
        dfq2 = -((q[0]*k)[:,1:] - (q[0]*k)[:,:-1])
        betha1 = (dfq1 - r13*dfq2)/(r11 - r13)
        betha3 = (-dfq1 + r11*dfq2)/(r11 - r13)
        dim = np.shape(q)
        dim = (dim[0], dim[1],dim[2]-1)
        W1 = np.zeros(dim)
        W3 = np.zeros(dim)
        W1[0] = betha1*r11
        W1[1] = betha1*1.
        W3[0] = betha3*r13
        W3[1] = betha3*1.
        Am = np.zeros(dim)
        Ap = np.zeros(dim)
        Am[:] = W1[:]
        Ap[:] = W3[:]
        return W1, W3, Am, Ap
    
    
    def Riemann_Elasticity_2Dy(q, rho, k, c):
        r11 = 1/np.sqrt(k*rho)[:-1,:]
        r13 = -1/np.sqrt(k*rho)[1:,:]
        dgq1 = -((q[2]/rho)[1:,:] - (q[2]/rho)[:-1,:])
        dgq3= -((q[0]*k)[1:,:] - (q[0]*k)[:-1,:])
        betha1 = (dgq1 - r13*dgq3)/(r11 - r13)
        betha3 = (-dgq1 + r11*dgq3)/(r11 - r13)
        dim = np.shape(q)
        dim = (dim[0], dim[1]-1,dim[2])
        W1 = np.zeros(dim)
        W3 = np.zeros(dim)
        W1[0] = betha1*r11
        W1[2] = betha1*1.
        W3[0] = betha3*r13
        W3[2] = betha3*1.
        Bm = np.zeros(dim)
        Bp = np.zeros(dim)
        Bm[:] = W1[:]
        Bp[:] = W3[:]
        return W1, W3, Bm, Bp
    
    def Transvese_Riemann_Elasticity_2Dx(q, rho, k, Ap, Am, c):
        r11 = (1/np.sqrt(k*rho))[:-2,1:]
        r13 = (-1/np.sqrt(k*rho))[1:-1:,1:]
        gamma1 = (Ap[0,1:-1,:] - r13*Ap[2,1:-1,:])/(r11 - r13)
        dim = np.shape(r11)
        BmAp = np.zeros((3,dim[0],dim[1]))
        BmAp[0] = -c[:-2,1:]*gamma1*r11
        BmAp[1] = -c[:-2,1:]*gamma1*0.
        BmAp[2] = -c[:-2,1:]*gamma1*1.
        
        r11 = (1/np.sqrt(k*rho))[1:-1,1:]
        r13 = (-1/np.sqrt(k*rho))[2:,1:]
        gamma3 = (-Ap[0,1:-1,:] + r11*Ap[2,1:-1,:])/(r11 - r13)
        dim = np.shape(r11)
        BpAp = np.zeros((3,dim[0],dim[1]))
        BpAp[0] = c[2:,1:]*gamma3*r13
        BpAp[1] = c[2:,1:]*gamma3*0.
        BpAp[2] = c[2:,1:]*gamma3*1.
        
        r11 = (1/np.sqrt(k*rho))[:-2,:-1]
        r13 = (-1/np.sqrt(k*rho))[1:-1,:-1]
        gamma1 = (Am[0,1:-1,:] - r13*Am[2,1:-1,:])/(r11 - r13)
        dim = np.shape(r11)
        BmAm = np.zeros((3,dim[0],dim[1]))
        BmAm[0] = -c[:-2,:-1]*gamma1*r11
        BmAm[1] = -c[:-2,:-1]*gamma1*0.
        BmAm[2] = -c[:-2,:-1]*gamma1*1.
        
        r11 = (1/np.sqrt(k*rho))[1:-1,:-1]
        r13 = (-1/np.sqrt(k*rho))[2:,:-1]
        gamma3 = (-Am[0,1:-1,:] + r11*Am[2,1:-1,:])/(r11 - r13)
        dim = np.shape(r11)
        BpAm = np.zeros((3,dim[0],dim[1]))
        BpAm[0] = c[2:,:-1]*gamma3*r13
        BpAm[1] = c[2:,:-1]*gamma3*0.
        BpAm[2] = c[2:,:-1]*gamma3*1.
        return BpAp, BmAp, BpAm, BmAm
    
    
    def Transvese_Riemann_Elasticity_2Dy(q, rho, k, Bp, Bm, c):
        r11 = (1/np.sqrt(k*rho))[1:,:-2]
        r13 = (-1/np.sqrt(k*rho))[1:,1:-1]
        gamma1 = (Bp[0,:,1:-1] - r13*Bp[1,:,1:-1])/(r11 - r13)
        dim = np.shape(r11)
        AmBp = np.zeros((3,dim[0],dim[1]))
        AmBp[0] = -c[1:,:-2]*gamma1*r11
        AmBp[1] = -c[1:,:-2]*gamma1*1.
        AmBp[2] = -c[1:,:-2]*gamma1*0.
        
        r11 = (1/np.sqrt(k*rho))[1:,1:-1]
        r13 = (-1/np.sqrt(k*rho))[1:,2:]
        gamma3 = (-Bp[0,:,1:-1] + r11*Bp[1,:,1:-1])/(r11 - r13)
        dim = np.shape(r11)
        ApBp = np.zeros((3,dim[0],dim[1]))
        ApBp[0] = c[1:,2:]*gamma3*r13
        ApBp[1] = c[1:,2:]*gamma3*1.
        ApBp[2] = c[1:,2:]*gamma3*0.
        
        r11 = (1/np.sqrt(k*rho))[:-1,:-2]
        r13 = (-1/np.sqrt(k*rho))[:-1,1:-1]
        gamma1 = (Bm[0,:,1:-1] - r13*Bm[1,:,1:-1])/(r11 - r13)
        dim = np.shape(r11)
        AmBm = np.zeros((3,dim[0],dim[1]))
        AmBm[0] = -c[:-1,:-2]*gamma1*r11
        AmBm[1] = -c[:-1,:-2]*gamma1*1.
        AmBm[2] = -c[:-1,:-2]*gamma1*0.
        
        r11 = (1/np.sqrt(k*rho))[:-1,1:-1]
        r13 = (-1/np.sqrt(k*rho))[:-1,2:]
        gamma3 = (-Bm[0,:,1:-1] + r11*Bm[1,:,1:-1])/(r11 - r13)
        dim = np.shape(r11)
        ApBm = np.zeros((3,dim[0],dim[1]))
        ApBm[0] = c[:-1,2:]*gamma3*r13
        ApBm[1] = c[:-1,2:]*gamma3*1.
        ApBm[2] = c[:-1,2:]*gamma3*0.
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
        F = np.zeros(np.shape(W1))
        norm = np.sum(np.square(W1[:,1:-1]),0)
        norm_zero_mask = np.array((norm == 0), dtype=float)
        norm_nonzero_mask = (1.0-norm_zero_mask)
    
        theta = np.ma.array(np.sum(W1[:,2:]*W1[:,1:-1],0))
        theta /= np.ma.array(norm)
        theta.fill_value = 0.
        theta = theta.filled()
        
        W1L = W1[:,1:-1]*norm_zero_mask + Limitador1D(theta, L)*W1[:,1:-1]*norm_nonzero_mask
        
        norm = np.sum(np.square(W2[:,1:-1]),0)
        norm_zero_mask = np.array((norm == 0), dtype=float)
        norm_nonzero_mask = (1.0-norm_zero_mask)
        
        theta = np.ma.array(np.sum(W2[:,0:-2]*W2[:,1:-1],0))
        theta /= np.ma.array(norm)
        theta.fill_value = 0.
        theta = theta.filled()
        
        W2L = W2[:,1:-1]*norm_zero_mask +  Limitador1D(theta, L)*W2[:,1:-1]*norm_nonzero_mask
        
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
            c = np.min([b , 2.*theta], 0)
            d = np.min([2.*b, theta], 0)
            phy = np.max([a,c,d],0)        
            return phy
    
        if o == 4:
            a = np.zeros(shape)
            b = np.ones(shape)
            c = np.min([(b + theta)/2. , 2.*theta, 2.*b], 0)
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
    
    Ibc = 3 #tipo de condicion de frontera a la izquierda
    Dbc = 3 #tipo de condicion de frontera a la derecha
    Arbc = 2 #tipo de condicion de frontera de la parte inferior del dominio 2D
    Abbc = 3 #tipo de condicion de frontera de la parte de arriba del dominio 2D
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
    seg = np.load('seg_eage_5409x1681.npy')
    # c = np.load('part_marmousi.npy')
    ref_factor1 = 3 - ref_factor
    c1 = seg[::2**ref_factor, ::2**ref_factor]
    
    Lim_Infx = 0; Lim_Supx = 13520; NumVolx = c1.shape[0]
    Lim_Infy = 0; Lim_Supy = 4200; NumVoly = c1.shape[1]
    t_inicial = 0.; t_final = 2000 ; n_salidas = t_final//200
    f0 = 0.010; # Frecuencia en Khrz
    CFL = 0.5 #Condicion CFL para el calculo del paso de tiempo
    Dimensional_Splitting = True # Si usamos division dimensional
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    Xc, Yc, Xn, Yn, dx, dy = Dominio(Lim_Infx, Lim_Supx, Lim_Infy, Lim_Supy, NumVolx, NumVoly) #Generar dominio computacional
    zsource = dy * (2**3//2**ref_factor) * 2 ; xsource = (Xc[0, Xc.shape[1]//2]) #Localizacion de la uente de ricker (debe ser un punto de la malla)
    yrec = zsource #profundidad del receptor 
    print (dy)
    print (dx)
    print (ref_factor)
    print (ref_factor1)
    #-----------------------------------------------------------------------------
    #Ricker-----------------------------------------------------------------------
    #Calculo de las velocidades del modelo (dominio)
    c = np.zeros((NumVolx + 4,NumVoly + 4))
    c[2:-2, 2:-2] = c1
    c = c.T
    c[:2,:] = c[2:3, :]
    c[-2:,:] = c[-4:-3,:]
    c[:,:2] = c[:,2:3]
    c[:,-2:] = c[:,-4:-3]
    
    #Calculo de la densidad y el modulo de compresibilidad de acuerdo a las velocidades
    b = 1 / (0.31 * (1e3*c)**0.25)
    b[c < 1.51] = 1.0
    Rho = 1./b
    K = c**2*Rho
    #-----------------------------------------------------------------------------
    VelMax = np.max(c) #Velocidad maxima para el calculo del CFL
    # dt = (np.minimum(dx,dy) * CFL)/VelMax #Tamaño del paso del tiempo
    dt = 0.5676
    q = Condicion_Inicial(Xc, Yc, Xn, Yn, NumVolx, NumVoly) #Calculo de la condicion inicial
    q = BCx(q, Ibc, 0)
    q = BCx(q, Dbc, 1)
    q = BCy(q, Arbc, 1)
    q = BCy(q, Abbc, 0)
    qb = np.copy(q) #Copia de q, para evolucion temporal
    Nt = int(round(t_final/dt))
    print(Nt)
    print(dt)
    receiver = np.zeros((Nt, NumVolx))
    inicio = time.time()
    if Dimensional_Splitting:
            for i in range(Nt):
                if i%(Nt//n_salidas) == 0:
                    print(t_inicial)
                    np.save("p_dswpa_t_"+str(int(t_inicial))+"_dx_" + str(int(dx)) + "_Test5.npy",(q[0]*K)[2:-2,2:-2])
                
                F = np.zeros((3,NumVoly+4,NumVolx+3))
                G = np.zeros((3,NumVoly+3,NumVolx+4))
                W1, W3, Am, Ap = Riemann_Elasticity_2Dx(q, Rho, K, c)
                qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (Ap[:,:,1:-2] + Am[:,:,2:-1])
                F[:,:,1:-1] = Fimm2D(c, W1, W3, Lim, 1)
                qb[:,:,2:-2] = qb[:,:,2:-2] - dt/dx * (F[:,:,2:-1] - F[:,:,1:-2])
                qb = BCx(qb, Ibc, 0)
                qb = BCx(qb, Dbc, 1)
                qb = BCy(qb, Arbc, 1)
                qb = BCy(qb, Abbc, 0)
                W1, W3, Bm, Bp = Riemann_Elasticity_2Dy(qb, Rho, K, c)
                qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (Bp[:,1:-2,:] + Bm[:,2:-1,:])
                G[:,1:-1,:] = Fimm2D(c, W1, W3, Lim, 2)
                qb[:,2:-2,:] = qb[:,2:-2,:] - dt/dy * (G[:,2:-1,:] - G[:,1:-2,:])
                
                t_inicial += dt
                receiver[i, :] = q[0,np.where(Yc[:,0]==yrec)[0][0],2:-2]
                qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2],f0, dt, t_inicial, Xc[2:-2,2:-2], Yc[2:-2, 2:-2], dx, dy, ref_factor1, yrec)
                qb = BCx(qb, Ibc, 0)
                qb = BCx(qb, Dbc, 1)
                qb = BCy(qb, Arbc, 1)
                qb = BCy(qb, Abbc, 0)
                
                q = np.copy(qb)
    else:
        for i in range(Nt):
            if i%(Nt//n_salidas) == 0:
                  np.save("p_fwpa_t"+str(int(t_inicial))+"_Test3.npy",(q[0]*K)[2:-2,2:-2])
                  print(t_inicial)
            F = np.zeros((3,NumVoly+4,NumVolx+3))
            G = np.zeros((3,NumVoly+3,NumVolx+4))
            W1x, W3x, Am, Ap = Riemann_Elasticity_2Dx(q, Rho, K, c)
            F[:,:,1:-1] = F[:,:,1:-1] + Fimm2D(c, W1x, W3x, Lim, 1)
            BpAp, BmAp, BpAm, BmAm = Transvese_Riemann_Elasticity_2Dx(q, Rho, K, Ap, Am, c)
            G[:,:-1,1:] = G[:,:-1,1:] - dt/(2.*dx)*BmAp
            G[:,1:,1:] = G[:,1:,1:] - dt/(2.*dx)*BpAp
            G[:,:-1,:-1] = G[:,:-1,:-1] - dt/(2.*dx)*BmAm
            G[:,1:,:-1] = G[:,1:,:-1] - dt/(2.*dx)*BpAm
            
            W1y, W3y, Bm, Bp = Riemann_Elasticity_2Dy(q, Rho, K, c)
            G[:,1:-1,:] = G[:,1:-1,:] + Fimm2D(c, W1y, W3y, Lim, 2)
            ApBp, AmBp, ApBm, AmBm = Transvese_Riemann_Elasticity_2Dy(q, Rho, K, Bp, Bm, c)
            F[:,1:,:-1] = F[:,1:,:-1] - dt/(2.*dy)*AmBp
            F[:,1:,1:] = F[:,1:,1:] - dt/(2.*dy)*ApBp
            F[:,:-1,:-1] = F[:,:-1,:-1] - dt/(2.*dy)*AmBm
            F[:,:-1,1:] = F[:,:-1,1:] - dt/(2.*dy)*ApBm
            
            
            qb[:,2:-2:,2:-2] = qb[:,2:-2,2:-2] - dt/dx * (Ap[:,2:-2,1:-2] + Am[:,2:-2,2:-1])\
                - dt/dy * (Bp[:,1:-2,2:-2] + Bm[:,2:-1,2:-2])\
                    - dt/dx * (F[:,2:-2,2:-1] - F[:,2:-2,1:-2])\
                        - dt/dy * (G[:,2:-1,2:-2] - G[:,1:-2,2:-2])
            
            
            t_inicial += dt
            receiver[i, :] = q[0,np.where(Yc[:,0]== yrec )[0][0],2:-2]
            qb[:,2:-2:,2:-2] = Ricker(qb[:,2:-2:,2:-2], f0, dt, t_inicial, Xc[2:-2,2:-2], Yc[2:-2, 2:-2], dx, dy, ref_factor1,yrec)                             
            
            
            qb = BCx(qb, Ibc, 0)
            qb = BCx(qb, Dbc, 1)
            qb = BCy(qb, Arbc, 1)
            qb = BCy(qb, Abbc, 0)
            #----------------------------------------------------------------------
            q = np.copy(qb)
    
    fin = time.time()
    print('Tiempo de ejecucion: ', fin - inicio)
    #plot = plt.imshow(c, cmap = 'jet', vmin = np.min(c), vmax = np.max(c))
    #ax = plt.gca()
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #cbar = plt.colorbar(plot, cax=cax)
    #cbar.set_label('Velocity (km/s)')
    
    # plt.figure(figsize=(15, 15))
    # plt.imshow((q[0]*K)[2:-2,2:-2], cmap = 'seismic', extent = [Lim_Infx, Lim_Supx, Lim_Supy, Lim_Infy])
    
    # plot1D_WPA = (q[0, :, NumVolx//2]*K[:, NumVolx//2])[2:-2]
    # yy_WPA = np.linspace(Lim_Infy, Lim_Supy, NumVoly)
    # fig, ax1 = plt.subplots()
    # cnt1 = ax1.plot(plot1D_WPA[::-1], yy_WPA[::-1])
    # plt.ylim(Lim_Supy, 0)
    
    #extent = [Lim_Infx, Lim_Supx, 1e-3*t_final, 0.]
    #aspect = Lim_Supx (1e-3*t_final)/.5
    
    #plt.figure(figsize=(15, 15))
    #plt.imshow(receiver, vmin=-0.01, vmax=0.01, cmap="seismic",
               # interpolation='lanczos', extent=extent, aspect=aspect)
    #plt.ylabel("Time (s)", fontsize=20)
    #plt.xlabel("Receiver position (m)", fontsize=20)
    if Dimensional_Splitting:
        np.save("p_dswpa_t_"+str(int(t_inicial))+"_dx_" + str(int(dx)) + "_Test5.npy",(q[0]*K)[2:-2,2:-2])
        np.save("rec_dswpa_t_2000_dx"+ str(int(dx)) +"_Test5.npy", receiver)
    else:
        np.save("p_fwpa_t"+str(int(t_inicial))+"_rf_" + str(ref_factor) + "_Test3.npy",(q[0]*K)[2:-2,2:-2])
        np.save("rec_fwpa_t1000_rf_"+ str(ref_factor) + "_Test3.npy", receiver)




ejecutar_dswpa(3)
ejecutar_dswpa(2)

