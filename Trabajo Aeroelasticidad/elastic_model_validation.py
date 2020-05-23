# -*- coding: utf-8 -*-
"""
******************* UNSTEADY AIRFOIL SOLVER    by: ENRIQUE FAUSTINO SILVELA PEÑA ***************************************
"""
#%% IMPORTAR LIBRERIAS? (PUES CLARO QUE ME IMPORTAN)
import numpy as np
import matplotlib.pyplot as plt
import Surface as Surface
from scipy.integrate import solve_ivp
from Grid import Grid



#%% UNIFORM FLOW
Vinf =1
alfa = 15
alfaR = alfa*(np.pi/180)
wo = 2.1

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.01
t = 0
nTime = 10000

#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)
airfoil = Surface.Line(0,-1,1,0,5)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% PARAMETROS ESTRUCTURALES
b = 1
ah = -0.5
xa = 0.2
ra = 0.5
wh = 0.2
wa = 0.2
mu = 100
Mass = mu*np.pi*rho*b*b

#%% ESTELA
Wake = []

#%% GRID
sizeX = [-2,2]
sizeY = [-1,2]
numX = 100
numY = 100
grid = Grid(sizeX,sizeY,numX,numY)

#%% KINEMATICA

# Velocidad del perfil, el origen se mueve en line recta na mas
theta = np.ones(nTime) *alfaR  
dtheta = np.zeros(nTime)

h = np.zeros(nTime)
dh = np.zeros(nTime)

dX0 = -np.ones(nTime)*Vinf
dY0 = -np.zeros(nTime)

#Unas transformaciones to flamas que salen en Plotkins
def RotationMatrix(alfa):
    return np.array([[np.cos(alfa), -np.sin(alfa)],
                    [np.sin(alfa),  np.cos(alfa)]])
    
    

#%% BEGIN TIME LOOP
Cl = np.zeros(nTime)
Cm = np.zeros(nTime)

Ut = np.zeros(nTime)
Wt = np.zeros(nTime)

#h[0] = 1
theta[0] = alfaR


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
   

#%% PLOT
#plt.cla()
#fig1 = plt.figure(0)
#PlotAirfoil(airfoil, Wake, grid, Ut[-1]- dtheta[-1]*h[-1],Wt[-1] - dh[-1] )
#PlotWake(Wake)

#%%

fig2 = plt.figure(2)
plt.title('theta0 = ' + str(alfa)+ 'º')
plt.xlabel('t')
plt.ylabel('Theta')


p1 = plt.plot(theta)
p2 = plt.plot(dtheta)
plt.legend((p1[0], p2[0]),('Theta','dTheta'))

fig3 = plt.figure(3)
plt.title('h0 = '+ str(h[0]))
plt.xlabel('t')
plt.ylabel('h')

p1 = plt.plot(h)
p2 = plt.plot(dh)
plt.legend((p1[0], p2[0]),('H','dH'))

