# -*- coding: utf-8 -*-
"""
******************* Validación del metodo no estacionario: Movimiento vertical ***************************************
 - Validación del código por medio de una aceleración instantánea, solución de Wagner
"""
#%% IMPORTAR LIBRERIAS? (PUES CLARO QUE ME IMPORTAN)
import numpy as np
import matplotlib.pyplot as plt
import Surface as Surface
from Grid import Grid
from ElementalSolutions import Vortex
from VORT2D import VORT2D



#%% UNIFORM FLOW
Vinf =1
alfa = 0
alfaR = alfa*(np.pi/180)
wo = 0

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.065
t = 0
nTime = 100

#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)
airfoil = Surface.Line(0,0,1,0,15)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% ESTELA
Wake = []
RemoveSmallVortices = 1

def UpdateWakePosition(Ut,Wt):
    for w in Wake:
        #Eliminamos vortices chiquitos?
        if(RemoveSmallVortices and (np.abs(w.gamma)<0.001 or w.x > 20)):
            Wake.remove(w)
            break
        vp, wp = 0,0 
        for p in airfoil.panels:
            a,b = VORT2D(p.gamma,w.x,w.y,p.xv,p.yv)
            vp += a
            wp += b
        vw, ww = 0, 0
        for w2 in Wake:
            if(w2!=w):
                vw = vw + w2.GetSpeedAt(w.x,w.y)[0]
                ww = ww + w2.GetSpeedAt(w.x,w.y)[1]
                
        w.x = w.x + (Ut + vw + vp)*dt
        w.y = w.y + (Wt + ww + wp)*dt

#%% GRID
sizeX = [-1,3]
sizeY = [-1,1]
numX = 150
numY = 150
grid = Grid(sizeX,sizeY,numX,numY)

#%% KINEMATICA

# Velocidad del perfil, el origen se mueve en line recta na mas
theta = np.ones(nTime) *alfaR  
dtheta = np.zeros(nTime)

h = np.zeros(nTime)
dh = -np.ones(nTime)*0.1
dh[0 :10] = 0

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

for i in range(nTime):
    #El tiempo se define aqui:
    t = dt*i

    #Integración simple de la velosidad       
    Ut[i],Wt[i] = RotationMatrix(theta[i]).dot([-dX0[i],-dY0[i]]) 

    # RESOLVEMOS LAS GAMMAS
    gamma = airfoil.GammaSolve(Ut[i],Wt[i],Wake,h[i],dh[i],dtheta[i],dt)    
    
    # RESOLVEMOS LA ESTELA
    # Primero la movemos:
    UpdateWakePosition(Ut[i],Wt[i])
    Wake.append(Vortex(gamma[-1],airfoil.trailX,0))
    
    # Calculamos fuerzas    
    L, M0 = airfoil.ForceSolve2(Ut[i],Wt[i],Wake)
    
    Cl[i] = L/(0.5*rho*(Vinf**2)*airfoil.c) *dt
    Cm[i] = M0/(0.5*rho*(Vinf**2)*(2**2))

#%% Solucion de kussner:
t = np.linspace(0,nTime*dt,nTime)*Vinf/(0.5)
Clest = 2*np.pi*(-dh/Vinf)
Clw = 2*np.pi*(-dh/Vinf)*(1-0.165*np.exp(-0.0455*t)-0.335*np.exp(-0.3*t))

fig1 = plt.figure(1)
plt.title('Cl - Wagner')
plt.xlabel('t')
plt.ylabel('Cl')


plt.plot(Cl,label = 'Cl numérico')
plt.plot(Clest,label = 'Cl estacionario')
plt.plot(Clw,label = 'Solución de Wagner')
plt.legend()