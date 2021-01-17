# -*- coding: utf-8 -*-
"""
******************* Validación del metodo no estacionario: Movimiento Oscilatorio ***************************************
 - 
"""
#%% IMPORTAR LIBRERIAS
import numpy as np
import matplotlib.pyplot as plt
import Surface as Surface
from Grid import Grid
from ElementalSolutions import Vortex
from VORT2D import VORT2D
from unsteady_plotter import PlotAirfoil, PlotWake



#%% UNIFORM FLOW
Vinf =1
alfa = 0
alfaR = alfa*(np.pi/180)
wo = 2

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.05
t = 0
nTime = 250

#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)
airfoil = Surface.Parabolic(0,0,1,0,15)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% ESTELA
Wake = []
RemoveSmallVortices = 0

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
sizeX = [-0.5,2]
sizeY = [-1,1]
numX = 50
numY = 50
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

for i in range(nTime):
    #El tiempo se define aqui:
    t = dt*i
    h[i] = 0.09*np.sin(wo*t)
    dh[i]= 0.09*wo*np.cos(wo*t)

    #Integración simple de la velosidad       
    Ut[i],Wt[i] = RotationMatrix(theta[i]).dot([-dX0[i],-dY0[i]]) 

    # RESOLVEMOS LAS GAMMAS
    gamma = airfoil.GammaSolve(Ut[i],Wt[i],Wake,h[i],dh[i],dtheta[i],dt)    
    
    # RESOLVEMOS LA ESTELA
    # Primero la movemos:
    UpdateWakePosition(Ut[i],Wt[i])
    Wake.append(Vortex(gamma[-1],airfoil.trailX,0))
    
    # Calculamos fuerzas
    L = Vinf*rho*sum([p.gamma for p in airfoil.panels]) 
    M0 = -rho*Vinf*sum([p.gamma*np.cos(p.phi)*p.xc for p in airfoil.panels])
    

    
    Cl[i] = L/(0.5*rho*(Vinf**2)*airfoil.c) 
    Cm[i] = M0/(0.5*rho*(Vinf**2)*(2**2))

   

#%% PLOT
plt.cla()
fig1 = plt.figure(0)
plt.title('Airfoil: Parabolic' + ', dt: ' + str(dt) + ', t: ' + str(int(dt*nTime)) + ', alfa0: ' + str(alfa) + ', wo: ' + str(wo) )
PlotAirfoil(airfoil, Wake, grid, Ut[-1]- dtheta[-1]*h[-1],Wt[-1] - dh[-1] )
PlotWake(Wake)

#%% Solucion de kussner:
t = np.linspace(0,nTime*dt,nTime)*Vinf/0.5
Clest = 2*np.pi*alfaR*np.ones(nTime)
Clk = 2*np.pi*alfaR*(1-0.5*np.exp(-0.130*t)-0.5*np.exp(-t))

fig1 = plt.figure(1)
plt.title('Cl')
plt.xlabel('t')
plt.ylabel('Cl')


plt.plot(Cl,label = 'Cl numérico')

plt.legend()

fig2 = plt.figure(2)
plt.title('Theta')
plt.xlabel('t')
plt.ylabel('Theta')


p1 = plt.plot(theta)
p2 = plt.plot(dtheta)
plt.legend((p1[0], p2[0]),('Theta','dTheta'))

fig3 = plt.figure(3)
plt.title('h')
plt.xlabel('t')
plt.ylabel('h')

p1 = plt.plot(h)
p2 = plt.plot(dh)
plt.legend((p1[0], p2[0]),('H','dH'))