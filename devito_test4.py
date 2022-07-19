#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                   as np
#import segyio
import matplotlib.pyplot       as plt
#==============================================================================

#==============================================================================
# Devito Imports
#==============================================================================
from devito import *
from examples.seismic.source import RickerSource, Receiver, TimeAxis
from examples.seismic.model import Model
#==============================================================================

#==============================================================================
# Latex Configuration
#==============================================================================
# init_printing(use_latex='mathjax')
# mpl.rcParams.update(mpl.rcParamsDefault)
#==============================================================================

#==============================================================================
# Plot Configuration
#==============================================================================
#==============================================================================

#==============================================================================
# plt.close("all")
#==============================================================================

#==============================================================================
def test3(ref_factor):
    
    rf         = ref_factor
    shape      = (int(1000 * 2**rf + 1), int(200*2**rf + 1))
    v          = np.empty(shape, dtype=np.float32)
    vel_media  = np.array([1.5, 1.81158421, 2.37165476, 3.4899143, 4.76228831, 5.231346, 6.3248646, 7., 7., 6.213656])
    n_layers   = 8
    layer_type = np.zeros(n_layers)
    prof       = np.array([0,-20,40,60,-60, -20,-5,10,10,10]) * 2**rf

    for k in range(n_layers):
        
        prof[k] = prof[k] + np.rint(shape[1]/n_layers)*k;

    i_vec = np.arange(0,shape[0])
    j_vec = np.arange(0,shape[1])

    eq_layer = np.zeros((n_layers, shape[0]))
    incl     = np.array([0.5,0.04,-0.07,-0.08,0.08,0.,0.,-0.07,0.,0.])

    for k in range(0,n_layers):
        
        if(layer_type[k]==0):
        
            eq_layer[k,:] = incl[k] * i_vec + prof[k]; 
              
    i_vec  = np.array([i_vec,]*shape[1]).transpose()
    j_vec  = np.array([j_vec,]*shape[0])
    v[:,:] = vel_media[0]

    for k in range(0, n_layers):

        indices = j_vec >= eq_layer[k, i_vec]
        
        v[np.where(indices)] = vel_media[k]
        
    return v
#==============================================================================

#==============================================================================
ref_factor = 1
#===========================================================================
#seg_vel = np.load("seg_eage_xcut_338.npy")

#marmousi_vel = np.load('part_marmousi.npy')
# c = np.load('part_marmousi.npy')
#c = marmousi_vel[::ref_factor, ::ref_factor]
#==============================================================================
# Setup Configuration
#==============================================================================
so        = 14
nptx      = int(100 * (2**ref_factor) + 1)
nptz      = int(100 * (2**ref_factor) + 1)
nptz1     = 0
x0        = 0.
x1        = 1000.
compx     = x1-x0
z0        = 0.   
z1        = 1000.
compz     = z1-z0
hxv       = compx/(nptx-1)
hzv       = compz/(nptz-1)
t0        = 0.0    # Milisencods
tn        = 250    # Milisencods
f0        = 0.015  # KHertz
nrec      = nptx
nbl       = 0      # int(0.5*nptx)
# so = 20
# nptx = 401
# nptz = 401
# nptz1 = 0
# x0 = 0
# x1 = 1000
# compx = (x1 - x0)
# hxv = (x1 - x0)/(nptx - 1)
# z0 = 0
# z1 = 1000
# compz = (z1 - z0)
# hzv = (z1 - z0)/(nptz - 1)
# t0 = 0
# tn = 250
# f0 = 0.015
# nrec = nptx
# nbl = 0

#==============================================================================

#==============================================================================
# Model Construction
#==============================================================================
c            = np.ones((nptx, nptz + nptz1))*1.5
# c[:, :nptz]  = test3(ref_factor)
# c[:, nptz:]  = c[:, nptz-1:nptz]

b            = 1 / (0.31 * (1e3*c)**0.25)
b[c < 1.51]  = 1.0

model        = Model(vp = c , b=b, shape=(nptx,nptz + nptz1),spacing=(hxv,hzv),nbl=nbl,space_order=so,origin=(x0,z0),extent=(compx,compz))
aspect_ratio = model.shape[0]/model.shape[1]
#==============================================================================

#==============================================================================
# Symbolic Dimensions
#==============================================================================
x,z  = model.grid.dimensions
zz_d = np.linspace(x0, x1, nptz)
t    = model.grid.stepping_dim
time = model.grid.time_dim
s    = time.spacing
#==============================================================================

#==============================================================================
#Parameters
#==============================================================================
ro  = 1./model.b
l2m = model.vp **2 * ro
#==============================================================================

#==============================================================================
# Time Construction
dt0        = model.critical_dt
time_range = TimeAxis(start=t0, stop=tn, step=dt0)
nt         = time_range.num - 1
# n_salidas  = int(tn/200) + 1
# jump       = nt//(n_salidas - 1)

#==============================================================================
# Ricker Source Construction
#==============================================================================
src                     = RickerSource(name='src',grid=model.grid,f0=f0,time_range=time_range)
xsource                 = 0.5 * (x1 - x0)
zsource                 = 0.5 * (z1 - z0)
src.coordinates.data[:] = [xsource,zsource]
#==============================================================================

#==============================================================================
# Symbolic Fields Construction
#==============================================================================
p = TimeFunction(name='p', grid = model.grid, staggered = NODE, space_order = so, time_order = 2)
v = VectorTimeFunction(name='v', grid = model.grid, space_order = so, time_order =2)
#==============================================================================

#==============================================================================
# Source Term Construction
#==============================================================================
src_p = src.inject(field = p.forward, expr = s * src * 4**(ref_factor))
#==============================================================================

#==============================================================================
# Receiver Term Construction
#==============================================================================
rec                       = Receiver(name="rec",grid=model.grid,npoint=nrec,time_range=time_range)
rec.coordinates.data[:,0] = np.linspace(x0,x1,num=nrec)
rec.coordinates.data[:,1] = 0.5 * (z1 - z0)
rec_term                  = rec.interpolate(expr=p)
#==============================================================================

#==============================================================================
# Symbolic Equation Construction
#==============================================================================
u_v = Eq(v.forward, v + s/ro * grad(p))
u_p = Eq(p.forward, p + s * l2m * div(v.forward))
#==============================================================================

#==============================================================================
# Save Plots for P
# time_subsampled = ConditionalDimension('t_sub',parent=time, factor=jump)
# psave = TimeFunction(name='psave',grid=model.grid,time_order=2,space_order=so,save=n_salidas,time_dim=time_subsampled,staggered=NODE)
#Pg    = np.zeros((n_salidas,nptx,nptz+nptz1))

#==============================================================================
# Operator Definition
# op2 = Operator([u_v, u_p] + src_p + rec_term + [Eq(psave,p)], subs=model.grid.spacing_map) 
op2 = Operator([u_v, u_p] + src_p + rec_term)

#==============================================================================
# Operator Evolution
#==============================================================================
op2(dt=dt0, time=nt)
# Pg[:]               = psave.data[:]
# Pg[n_salidas-1,:,:] = p.data[0,:,:]
#==============================================================================

  # Graphical Plots
  # plt.figure()
  
# plt.figure()
plt.imshow(p.data[0,:,:].T, cmap = 'seismic', extent = [x0, x1, z1, z0])
plt.show()

# # Graphical Plots for Save Steptimes
# plt.figure()
# plt.imshow(rec.data, cmap = 'seismic')
# plt.show()
#==============================================================================
# for i in range(1,n_salidas):
    # np.save('p_devito_ref_dx' + str(int(hxv)) +'_SO_'+ str(so) +'_t_' + str(int(i*200)) +'_Test4.npy',psave.data[i])
    # fscale = 10**(-3)
    # extent = [fscale*x0,fscale*x1, fscale*z1, fscale*z0]
    # plt.title('Psave = %d'%(i))
    # plt.imshow(np.transpose(Pg[i,:,:]), cmap = 'seismic', extent=extent)
    # plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    # plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    # plt.savefig('psave_%d.png'%(i),dpi=100)
    # plt.show()
    # plt.close()
#==============================================================================
# np.save('rec_devito_dx_'+str(int(hxv))+'_SO_'+ str(so) +'_t_' + str(tn) + '_Test4.npy',rec.data)
#==============================================================================
#np.save('p_devito_referencia_Test3.npy',p.data[0])
#np.save('rec_devito_referencia_Test3.npy',rec.data)
# ==============================================================================
