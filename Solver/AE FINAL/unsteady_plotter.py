import numpy as np
import matplotlib.pyplot as plt
from VORT2D import VORT2D


def PlotAirfoil(airfoil,Wake,grid, dX,dY):
    plt.gca().set_aspect('equal')   
    Cp = np.zeros([grid.numX,grid.numY])
    Vx, Vy = np.zeros([grid.numX,grid.numY]), np.zeros([grid.numX,grid.numY])
    for i in range(grid.numX):
        for j in range(grid.numY):
            nu = 0
            nv = 0
            uw = 0
            vw = 0
            for p in airfoil.panels:
                u,v = VORT2D(p.gamma,grid.XX[i,j],grid.YY[i,j],p.xv,p.yv)
                nu = u + nu 
                nv = v + nv 
            for w in Wake:
                uw += w.GetSpeedAt(grid.XX[i,j],grid.YY[i,j])[0]
                vw += w.GetSpeedAt(grid.XX[i,j],grid.YY[i,j])[1]
            Vx[i,j] = Vx[i,j] + nu + uw
            Vy[i,j] = Vy[i,j] + nv + vw
    Vx = Vx + dX
    Vy = Vy + dY
    
    Vt = np.sqrt(Vx**2+Vy**2)
    Cp = 1- (Vt/dX)**2
#    plt.contourf(grid.XX,grid.YY,Cp,500,cmap='jet') 
    airfoil.PlotSurface()
#    grid.PlotStream(Vx,Vy)
    #airfoil.DrawShape()
    
def PlotWake(Wake):
    for i,w in enumerate(Wake):
        if(i<len(Wake)-1):
            plt.plot([w.x,Wake[i+1].x],[w.y,Wake[i+1].y],'k',linewidth = 0.5)  