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

#%%
Γ = np.array([[1,2],[2,1]]) 
R = np.sqrt(Γ)
# print(R)
# print(np.dot(R,R))

Γ = np.array([[2,1],[1,2]])
s = np.sqrt(Γ[0][0]*Γ[1][1]-Γ[1][0]**2)
t = np.sqrt(Γ[0][0]+Γ[1][1]+2*s)
sqrt_Γ = (1/t)*(Γ+s*np.identity(2))
print(Γ)
print(sqrt_Γ)
print(np.dot(sqrt_Γ,sqrt_Γ))
#%%

Γ = np.array([[1,2],[3,4]]) 
print(Γ[0][1])
# print(Γ[0])

#%% Funciones que probablemente no use

def lineal(x, a, b):
    return a*x+b

def gaussiana(r, sig):
#    r = np.sqrt((x-x0)**2 + (y-y0)**2)
    return np.e**(-(r**2)/(2*sig**2))/(sig*np.sqrt(2*np.pi))

def theta(x,y):
    r = np.sqrt(x**2 + y**2)
    if r == 0:
        theta = 0
    elif y >= 0 and r != 0:
        theta = np.arccos(x/r)
    else:
        theta = -1*np.arccos(x/r)
    return theta

def gaussian_force(x, y, x0, y0, sig, f_max): # es el resultante de un potencial gaussiano *(-1)
    r = np.sqrt((x-x0)**2 + (y-y0)**2)        # sig = std del potencial
    tita = theta(x-x0, y-y0)                  # máx valor de la fuerza 
    F_r = -f_max*(np.exp(0.5))*r*(np.exp(-(r**2)/(2*sig**2)))/sig
    F_x = F_r*np.cos(tita)
    F_y = F_r*np.sin(tita)
    return F_x, F_y

#%% Codigo de Simulación 
m = 0.2 
beta = 100

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


x0, y0 = 0, 0
vx_0, vy_0 = 5, 0
T=1
N = 1000

x1, y1 = simu(x0,y0,vx_0,vy_0,T,N)
# print(x1)


plt.close("all")
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.scatter(x1[0], y1[0], color = "darkgreen", label = "Start", marker = "X",s = 100, zorder = 2)
ax.scatter(x1[-1], y1[-1], color = "darkblue", label = "End", marker = "X",s = 100, zorder = 2)
ax.axis('equal')


s = 3 #int(n/T_frames)
for i in range(N+1):           
#    ax.plot(x5[i:i+s+1], y5[i:i+s+1], linewidths=0, marker='o', s=3, cmap=plt.cm.winter, zorder = 1) # la traza    
    ax.plot(x1[i:i+2], y1[i:i+2], color="k", alpha= (i/N)) # la traza 
# ax.set_axis_off()
ax.legend()
plt.show()













































