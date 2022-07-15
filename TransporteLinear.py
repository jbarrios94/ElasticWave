import numpy as np
import matplotlib.pyplot as plt
def Exacta(x,n,t):
    u0 = np.zeros(n+1)
    for i in range(n+1):
		if (0.1 + t < x[i] < 0.5 + t):
			u0[i] = 1.
    return u0

def Dominio(inf,sup,n,dx):
    x = np.zeros(n+1)
    x[0] = inf
    for i in range(1,n+1):
        x[i] = x[i-1] + dx
    return x
    
def Fluxo(op,u):
	if op == 0:
		v = 1.0
		return v * u
	elif op == 1:
		
inf = 0.
N = 100
Cr = 0.2
dx = (sup-inf)/N
dt = (dx*Cr)/(v*MaxDer)
Nt = int(round(T/dt))
print(Nt)
x = Dominio(inf,sup,N,dx)
#u0 = Exacta(x,N,0)
u0 = np.zeros(N+1)
"""
Lax-Friedrichs
"""

u = np.copy(u0)
for i in range(Nt):
	term1 = 1./2.*(u0[0:N-2] + u0[2:N])
	term2 = 1./2.*dt/dx*(Fluxo(2,u0[2:N]) - Fluxo(2,u0[0:N-2]))
	u[1:N-1] = term1 - term2
	u[0] = 1
	u[N] = u0[N-1]
	u0 = np.copy(u)
	
plt.plot(x,u0)
plt.show()
return u*u/2.
	elif op == 2:
		a = 1./2.
		v = 1.
		phi = 1.
		return u*u/(u*u + a*(1.-u)*(1-u))
		
		
	
v = 1.
MaxDer = 2.5
T = 0.5
sup = 1.
inf = 0.
N = 100
Cr = 0.2
dx = (sup-inf)/N
dt = (dx*Cr)/(v*MaxDer)
Nt = int(round(T/dt))
print(Nt)
x = Dominio(inf,sup,N,dx)
#u0 = Exacta(x,N,0)
u0 = np.zeros(N+1)
"""
Lax-Friedrichs
"""

u = np.copy(u0)
for i in range(Nt):
	term1 = 1./2.*(u0[0:N-2] + u0[2:N])
	term2 = 1./2.*dt/dx*(Fluxo(2,u0[2:N]) - Fluxo(2,u0[0:N-2]))
	u[1:N-1] = term1 - term2
	u[0] = 1
	u[N] = u0[N-1]
	u0 = np.copy(u)
	
plt.plot(x,u0)
plt.show()
