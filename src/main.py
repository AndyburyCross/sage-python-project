from sage.all import *
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

# Sage symbolic maths
x=var('x')
print("Integral of sin(x)^2:")
print(integrate(sin(x)**2,x))

#Numerical ODE
def f(t, y):
	return -y + sin(t)

sol = solve_ivp(f, (0,10), [1.0])

plt.plot(sol.t, sol.y[0])
plt.xlabel("t")
plt.ylabel("y")
plt.title("ODE solution")
plt.savefig("output.png")

print("Saved plot to output.png")
