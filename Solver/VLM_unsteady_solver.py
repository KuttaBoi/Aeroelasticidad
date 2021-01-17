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
from ElementalSolutions import Vortex
from VORT2D import VORT2D
from unsteady_plotter import PlotAirfoil, PlotWake



#%% UNIFORM FLOW
Vinf =1
alfa = 0
alfaR = alfa*(np.pi/180)
wo = 2.1

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.05
t = 0
nTime = 500

#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)
airfoil = Surface.Line(0,0,2,0,15)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% PARAMETROS ESTRUCTURALES
b = 1
ah = -0.5
xa = 0.25
ra = 0.5
wh = 1.8
wa = 1
mu = 100
Mass = mu*np.pi*rho*b*b

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
            n,m = VORT2D(p.gamma,w.x,w.y,p.xv,p.yv)
            vp += n
            wp += m
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

h = np.zeros(nTime)
h[0] = 1
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



#%% Loopeamos cada instante de tiempo

def Lagrange(t,y,h,alfa,L,M0):
    A = [[1/b ,xa   ],
         [xa,b*ra**2]]
    
    bs = [[-L/(Mass*b) - h*(wh**2)/b],
          [M0/(Mass*b**2) - alfa*(ra*wa**2)/b**2]]
    return np.linalg.solve(A,bs)


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
#    L = Ut[i]*rho*sum([p.gamma for p in airfoil.panels]) 
#    M0 = rho*Ut[i]*sum([p.gamma*np.cos(p.phi)*p.xc for p in airfoil.panels])
    
    L, M0 = airfoil.ForceSolve(Ut[i],Wt[i],dtheta[i],dh[i],h[i],Wake)

    Cl[i] = L/(0.5*rho*(Vinf**2)*airfoil.c) 
    Cm[i] = M0/(0.5*rho*(Vinf**2)*(2**2))
#    
    # Solve H
    y0 = [dh[i],dtheta[i]]

    y = solve_ivp(Lagrange,[t,t+dt],y0,method = 'RK45',t_eval = [t+dt],vectorized = True,args=(h[i],theta[i],L,M0))
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
fig1 = plt.figure(1)
plt.title('Cl')
plt.xlabel('t')
plt.ylabel('Cl')


p1 = plt.plot(Cl,label = 'Cl')
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
