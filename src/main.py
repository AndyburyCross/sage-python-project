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

show(V(M,phi))
