# -*- coding: utf-8 -*-
"""
******************* Respuesta en Flutter ***************************************
"""
#%% IMPORTAR LIBRERIAS? (PUES CLARO QUE ME IMPORTAN)
import numpy as np
import matplotlib.pyplot as plt
import Surface as Surface
from Grid import Grid
from ElementalSolutions import Vortex
from VORT2D import VORT2D
from unsteady_plotter import PlotAirfoil, PlotWake


#%% UNIFORM FLOW
Vinf = 9
alfa = 5
alfaR = alfa*(np.pi/180)

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.1
t = 0
nTime = 1000

#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)
airfoil = Surface.Line(0,0,2,0,5)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% PARAMETROS ESTRUCTURALES
b = 1
ah = -0.5
xa = 0.25
ra = 0.5
wh = 0.2
wa = 1
mu = 100
Mass = mu*np.pi*rho*b*b

#%% ESTELA
Wake = []
RemoveSmallVortices = 1
maxWakeVortex = 25

def UpdateWakePosition(Ut,Wt):
    for w in Wake:
        #Eliminamos vortices chiquitos?
        if(RemoveSmallVortices and np.abs(w.gamma)<0.001):
            Wake.remove(w)
            break
        for p in airfoil.panels:
            vp,wp = VORT2D(p.gamma,w.x,w.y,p.xv,p.yv)
        vw, ww = 0, 0
        for w2 in Wake:
            if(w2!=w):
                vw = vw + w2.GetSpeedAt(w.x,w.y)[0]
                ww = ww + w2.GetSpeedAt(w.x,w.y)[1]
                
        w.x = w.x + (Ut + vw + vp)*dt
        w.y = w.y + (Wt + ww + wp)*dt

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
dda = np.zeros(nTime)

h = np.zeros(nTime)
dh = np.zeros(nTime)
ddh = np.zeros(nTime)

dX0 = -np.ones(nTime)*-Vinf*np.cos(alfaR)
dY0 = -np.zeros(nTime)*-Vinf*np.sin(alfaR)

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
    
    h[i] = dh[i]*dt + h[i]
    theta[i] = dtheta[i]*dt+ theta[i] -dh[i]/Vinf
    #Integración simple de la velosidad       
    Ut[i],Wt[i] = RotationMatrix(theta[i]).dot([-dX0[i],-dY0[i]]) 
    Ut[i] = Ut[i]
    Wt[i] = Wt[i] - dh[i] + dtheta[i]*(xa+ah)*b
    # RESOLVEMOS LAS GAMMAS
    gamma = airfoil.GammaSolve2(Ut[i],Wt[i],Wake)     
    
    # RESOLVEMOS LA ESTELA
    # Primero la movemos:
    UpdateWakePosition(Ut[i],Wt[i])
    if(len(Wake)>maxWakeVortex):
        Wake.pop(0)
    Wake.append(Vortex(gamma[-1],2.1,0))
    
    # Calculamos fuerzas
    L = rho*Vinf*sum([p.gamma for p in airfoil.panels]) + rho*np.pi*(b**2)*(ddh[i] + Vinf*dtheta[i])
    M0 = rho*Vinf*sum([p.gamma*(b*0.5-p.xc)*np.cos(theta[i]) for p in airfoil.panels])
    
    Cl[i] = L/(rho*(Vinf**2)*b) 
    Cm[i] = M0/(0.5*rho*(Vinf**2)*((b/2)**2))
    
    A = [[1/b,xa],[xa,(ra**2)/b]]
    bs = [[-L/(b*Mass) - (wh**2)*h[i]],[M0/(b*b*Mass) - theta[i]*((ra*wa)**2)]]
    
    if(i<nTime-1):
        ddh[i],dda[i] = np.linalg.solve(A,bs)
        dh[i+1] = ddh[i]*dt + dh[i]
        dtheta[i+1] = dda[i]*dt + dtheta[i]


#%% PLOT
#plt.cla()
#fig1 = plt.figure(0)
#PlotAirfoil(airfoil, Wake, grid, Ut[-1]- dtheta[-1]*h[-1],Wt[-1] - dh[-1] )
#PlotWake(Wake)

#%%
fig1 = plt.figure(1)
plt.title('Cl- Vinf: ' + str(Vinf))
plt.xlabel('t')
plt.ylabel('Cl')


p1 = plt.plot(Cl,label = 'Cl' )
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



