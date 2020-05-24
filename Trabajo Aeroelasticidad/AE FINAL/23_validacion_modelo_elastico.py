# -*- coding: utf-8 -*-
"""
******************* Validacion del método elástico ***************************************
"""
#%% IMPORTAR LIBRERIAS? (PUES CLARO QUE ME IMPORTAN)
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#%% INTEGRACION TEMPORAL
dt = 0.01
t = 0
nTime = 10000

#%% PARAMETROS ESTRUCTURALES
b = 1
ah = -0.5
xa = 0.2
ra = 0.5
wh = 0.2
wa = 0.2
mu = 100
Mass = mu*np.pi*1.22*b*b


#%% KINEMATICA

theta = np.ones(nTime)
dtheta = np.zeros(nTime)

h = np.zeros(nTime)
dh = np.zeros(nTime)

h[0] = 1
theta[0] = 0


#%% Loopeamos cada instante de tiempo

def Lagrange(t,y,h,alfa,L,M0):
    A = [[1 ,xa   ],
         [xa,ra**2]]
    
    bs = [[-L/Mass - h*(wh**2)],
          [M0/Mass - alfa*(wa**2)]]
    return np.linalg.solve(A,bs)


for i in range(nTime):
    #El tiempo se define aqui:
    t = dt*i
    y0 = [dh[i],dtheta[i]]

    y = solve_ivp(Lagrange,[t,t+dt],y0,method = 'RK45',t_eval = [t+dt],vectorized = True,args=(h[i],theta[i],0,0))
    if(i<nTime-1):
        dh[i+1] = y.y[0]
        dtheta[i+1] = y.y[1]
        h[i+1] = dh[i]*dt + h[i]
        theta[i+1] = dtheta[i]*dt+ theta[i]


#%%

fig2 = plt.figure(2)
plt.title('theta0 = ' + str(theta[0])+ 'º')
plt.xlabel('t')
plt.ylabel('Theta')

plt.plot(theta,label = 'theta')
plt.plot(dtheta, label = 'dtheta')
plt.legend()

fig3 = plt.figure(3)
plt.title('h0 = '+ str(h[0]))
plt.xlabel('t')
plt.ylabel('h')

plt.plot(h, label = 'h')
plt.plot(dh, label = 'dh')
plt.legend()


