#%% Librerias
"""
@author: Tucu

Copyright © 2024 Tucu. All rights reserved.

===============================
elTucu-labs
===============================

This product includes software developed by elTucu-labs. Project for use in the thesis project.

Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met:

*Redistributions of source code must remain the above copyright notice, this list of conditions and the
following disclaimer. 

*Redistributions in binary form must reproduce the above copyright notice, this list of conditions and 
the following disclaimer in the documentation and/or other materials provided with the distribution.


*Neither the name of the Lab nor the names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

"""

import numpy as np
import matplotlib.pyplot as plt
import random
import imageio as imageio            
from tqdm import tqdm
from IPython import get_ipython
get_ipython().run_line_magic("matplotlib","qt5")


#%% Codigo de Simulación 
m = 0.2 
beta = 10

def matrix_sqrt(A):
    s = np.sqrt(A[0][0]*A[1][1]-A[1][0]**2)
    t = np.sqrt(A[0][0]+A[1][1]+2*s)
    sqrt_A = (1/t)*(A+s*np.identity(2))
    return sqrt_A

def Mu_0(x,y): # = (f(x,y) - ∂xΦ(x,y) - ∂yΦ(x,y))/m
    m1 = 1 
    m2 = 1.5
    m3 = 1 
    fx = - m3*y
    fy = m3*x
    mu_0 = np.array([-(4*m1*x**3-2*m2*x)+fx,-2*y+fy])/m
    # mu_0 = [0,0]
    return mu_0 # = [mu_0x, mu_0y]

def Gamma(x,y):
    Γ = [[0.5,0.2],[0.2,0.5]]
    # Γ = [[1,0],[0,1]]
    # Γ = [[0.5,0],[0,0.5]]
    return Γ

def simu(x0,y0,vx_0,vy_0,T,N):
    dt = T/N
    X = [x0]
    Y = [y0]
    vx = vx_0
    vy = vy_0
    a1 = np.sqrt(2/(m*beta))
    for i in range(N-1):
        mu_0 = Mu_0(X[-1],Y[-1])
        gamma = Gamma(X[-1],Y[-1])
        sqrt_Γ = matrix_sqrt(gamma)
        eta_x , eta_y = np.random.normal(0, 1, 2) # mean, sigma, results
        X.append(X[-1] + vx*dt)
        Y.append(Y[-1] + vy*dt)
        vx = mu_0[0] - gamma[0][0]*vx - gamma[0][1]*vy + a1*sqrt_Γ[0][0]*eta_x + a1*sqrt_Γ[0][1]*eta_y
        vy = mu_0[1] - gamma[1][0]*vx - gamma[1][1]*vy + a1*sqrt_Γ[1][0]*eta_x + a1*sqrt_Γ[1][1]*eta_y 
    return X, Y

# Parametros para la simulacion

x0, y0 = 0, 0
vx_0, vy_0 = 5, 0
T=2 
N = 1000

x1, y1 = simu(x0,y0,vx_0,vy_0,T,N)
# print(x1)


plt.close("all")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(x1[0], y1[0], color = "darkgreen", label = "Start", marker = "X",s = 100, zorder = 2)
ax.scatter(x1[-1], y1[-1], color = "darkblue", label = "End", marker = "X",s = 100, zorder = 2)
ax.axis('equal')



for i in range(N+1):           
    ax.plot(x1[i:i+2], y1[i:i+2], color="k", alpha= (i/N)) # la traza 
# ax.set_axis_off()
ax.legend()
plt.show()

#%% Campo de fuerza: Vectores

# Meshgrid 
x, y = np.meshgrid(np.linspace(-1.5, 1.5, 10),  
                   np.linspace(-1, 1, 10)) 
  
# Directional vectors 
u = Mu_0(x, y)[0]*m
v = Mu_0(x, y)[1]*m
  

# Plotting Vector Field with QUIVER 
plt.close("all")
plt.quiver(x, y, u, v, color='darkblue') 
plt.title('Vector Field') 
  
# Setting x, y boundary limits 
plt.xlim(-1.5, 1.5) 
plt.ylim(-1, 1) 
  
# Show plot with grid 
plt.grid() 
plt.show() 


#%% Campo de fuerza: Lineas de fuerza


plt.close("all")

# Particula con traza
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(x1[0], y1[0], color = "darkgreen", label = "Start", marker = "X",s = 300, zorder = 4)
ax.scatter(x1[-1], y1[-1], color = "darkblue", label = "End", marker = "X",s = 300, zorder = 4)
ax.axis('equal')
for i in range(N+1):           
    ax.plot(x1[i:i+2], y1[i:i+2], color="k", alpha= (i/N),linewidth=2) # la traza 

# Campo de fuerza
# 1D arrays 
x = np.arange(-1,1,0.1) 
y = np.arange(-1,1,0.1) 
# Meshgrid 
X,Y = np.meshgrid(x,y) 
# Assign vector directions 
u = Mu_0(X, Y)[0]*m
v = Mu_0(X, Y)[1]*m
plt.streamplot(X,Y,u,v, density=1, linewidth=None, color='#A23BEC', arrowsize=2) 
  
# Show plot with grid 
plt.legend()
plt.grid() 
plt.show()








































