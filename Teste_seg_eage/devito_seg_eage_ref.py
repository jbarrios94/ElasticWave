#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy as np
#==============================================================================

#==============================================================================
# Devito Imports
from devito import *
from examples.seismic.source import RickerSource, Receiver, TimeAxis
from examples.seismic import plot_image, demo_model, Model
import Graficar as gf
import matplotlib.pyplot as plt

#==============================================================================
# Latex Configuration
#==============================================================================

#==============================================================================
# Setup Configuration
#==============================================================================
# c = np.load("/hom:e/juan/ElasticWave/Teste_seg_eage/seg_teste5_dx20.npy")
# c = c[::2**3, ::2**3]

# print(c.shape)
so = 8
# nptx = c.shape[0]
# nptz = c.shape[1]
nptx = 350 * 2**0 + 1
nptz = 210 * 2**0 + 1
c = np.ones((nptx, nptz)) * 1.5
x0 = 3000
x1 = 10000
compx = (x1 - x0)
hxv = (x1 - x0)/(nptx - 1)
print(hxv)
z0 = 0
z1 = 4200
compz = (z1 - z0)
hzv = (z1 - z0)/(nptz - 1)
print(hzv)
t0 = 0
tn = 1500
f0 = 0.005
nrec = nptx
nbl = 0
#==============================================================================

# =============================================================================
# Model Construction
# =============================================================================
b            = 1 / (0.31 * (1e3*c)**0.25)
b[c < 1.51]  = 1.0
model        = Model(vp = c , b=b, shape=(nptx,nptz),spacing=(hxv,hzv),nbl=nbl,space_order=so,origin=(x0,z0),extent=(compx,compz))
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
#==============================================================================
dt0        = model.critical_dt
print(dt0)
time_range = TimeAxis(start=t0, stop=tn, step=dt0)
nt         = time_range.num - 1
n_salidas  = int(tn/200) + 1
jump       = nt//(n_salidas - 1)
#==============================================================================

#==============================================================================
# Ricker Source Construction
#==============================================================================
src                     = RickerSource(name='src',grid=model.grid,f0=f0,time_range=time_range)
xsource                 = 0.5 * (x1 + x0)
zsource                 = 1400
src.coordinates.data[:] = [xsource,zsource]
# #==============================================================================

#==============================================================================
# Symbolic Fields Construction
#==============================================================================
p = TimeFunction(name='p', grid = model.grid, staggered = NODE, space_order = so, time_order = 2)
v = VectorTimeFunction(name='v', grid = model.grid, space_order = so, time_order =2)
#==============================================================================

#==============================================================================
# Source Term Construction
#==============================================================================
src_p = src.inject(field = p.forward, expr = s * src/ (hxv * hzv))
#==============================================================================

#==============================================================================
# Receiver Term Construction
#==============================================================================
rec                       = Receiver(name="rec",grid=model.grid,npoint=nrec,time_range=time_range)
rec.coordinates.data[:,0] = np.linspace(x0,x1,num=nrec)
rec.coordinates.data[:,1] = 1400
rec_term                  = rec.interpolate(expr=p)
#==============================================================================

#==============================================================================
# Symbolic Equation Construction
#==============================================================================
u_p = Eq(p.forward, p - s * l2m * div(v))
u_v = Eq(v.forward, v - s/ro * grad(p.forward))
#==============================================================================

#==============================================================================
# Save Plots for P
time_subsampled = ConditionalDimension('t_sub',parent=time, factor=jump)
psave = TimeFunction(name='psave',grid=model.grid,time_order=2,space_order=so,save=n_salidas,time_dim=time_subsampled,staggered=NODE)
Pg    = np.zeros((n_salidas,nptx,nptz))

#==============================================================================
# Operator Definition
op2 = Operator([u_p, u_v] + src_p + rec_term + [Eq(psave,p)], subs=model.grid.spacing_map) 
# op2 = Operator([u_v, u_p] + src_p + rec_term)

#==============================================================================
# Operator Evolution
#==============================================================================
op2(dt=dt0, time=nt)
#==============================================================================

# ==============================================================================
# for i in range(1,n_salidas):
#     np.save('p_devito_so8_t' + str(int(i*200)) +'_Test5.npy',psave.data[i])
# np.save('rec_devito_so8_t' + str(tn) + '_Test5.npy',rec.data)
# ==============================================================================
gf.plot_image(p.data[0].T, extent = [x0, x1, z1, z0])
gf.plot_seismic_traces([p.data[0, 180, :]], 0, 1500)
print(np.max(p.data[0, 180, :]))
plt.show()