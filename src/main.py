from sage.all import *
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import numpy as np

# a) Define variables used in Inverse Power Law:

phi = var(r'φ') # scalar field
k = 1
#Mpl = 1 # Reduced planck mass units
alp = var(r'α') # steepness of inverse power law
alp = 4
M = var('M') # Mass scale of IPL
N = var('N') # N units
beta = var(r'β') # coupling variable

# Expression for V and derivative
def V(M,phi):
    return (M**(4+alp))/phi**alp

def dVdphi(M,phi):
    return -alp*V(M,phi)/phi

#show(V(M,phi))

h = 0.7
ohm_m = 0.31 # Current approximate value for Ωm 
Ho = 2.1332e-42*h*4.10677e-19 # Hubble constant in reduced plank mass units
rho_c0 = (3*Ho**2/(k**2)) # current critical energy density
M_val = ((rho_c0*(8*np.pi)**(alp/2))**(1/(alp+4)))*.65 # Value set for M
show(M_val)
gam = 1
beta = 0
Ni = 0 

def CoupledSteinhardt(X,N):
    a = exp(N) # N = ln(a)
    ai = exp(Ni)
    rho_m = rho_c0*ohm_m*(a*ai)**(-3)
    Vfunc = V(M_val,X[0])
    #Vfunc = V.subs(V0 == V_value,phi ==X[0])
    dVfunc = dVdphi(M_val,X[0])

    H2 = (rho_m + Vfunc)/(3-0.5*(X[1])**2)
    HdotH2 = -0.5*((X[1])**2+ gam*rho_m/H2)
    x1prime = X[1]
    x2prime = -(3+HdotH2)*X[1] - dVfunc/H2- beta*rho_m/H2
    #show(x2prime)
    return [x1prime,x2prime]

Ni=0
show(CoupledSteinhardt([1,0],0))

Nrange = np.arange(-5,5,.1)
Ni =0
sol = odeint(CoupledSteinhardt,[1,0],t=Nrange)
# Plot of φ against redshift
phi_sol = sol[:,0]
phi_prime = sol[:,1]
#Vnumber = sol[:,2]
#show(Vnumber)
z_values = exp(-Nrange)
plt.xlabel("N")
plt.ylabel(r"φ $[M_{pl}]$")
plt.plot(Nrange,phi_sol,'--')
plt.savefig("phi.png")
