import numpy as np
from VORT2D import VORT2D
class Vortex:
    x = 0
    y = 0
    gamma = 1

    def __init__(self,gamma,x,y):
        self.x = x
        self.y = y
        self.gamma = gamma
        
    def GetSpeedAt(self,x,y):
        return VORT2D(self.gamma,x,y,self.x,self.y)
    
        

class Uniform:
    Vx = 0
    Vy = 0
    V  = 0
    Vinf = 0
    alpha = 0
    alphaR = 0
    def __init__(self,Vinf,alpha,XX,YY):
        sizeX = np.size(XX,1);
        sizeY = np.size(YY,1);
        self.Vx = np.ones([sizeX,sizeY]);
        self.Vy = np.ones([sizeX,sizeY]);
        self.V = np.ones([sizeX,sizeY]);
        
        self.Vx = self.Vx*Vinf*np.cos(alpha*(np.pi/180))
        self.Vy = self.Vy*Vinf*np.sin(alpha*(np.pi/180))
        self.V = np.sqrt(self.Vx**2 + self.Vy**2)
        self.Vinf = Vinf
        self.alpha = alpha*(np.pi/180)
        
        
    def GetPotentialAtGrid(self,XX,YY):
        numX = np.size(XX,1);
        numY = np.size(YY,1);
        
        Phi = np.zeros([numX,numY])
        for i in range(numX):
            for j in range(numY):
                Phi[i,j] = self.Vinf*np.cos(self.alpha)+self.Vinf*np.sin(self.alpha)
        return Phi    

        
 
     

class Source:
    Vx = 0
    Vy = 0
    V  = 0
    def __init__(self,lmbda,x0,y0,XX,YY):
        sizeX = np.size(XX,1);
        sizeY = np.size(YY,1);
        self.Vx = np.zeros([sizeX,sizeY]);
        self.Vy = np.zeros([sizeX,sizeY]);
        self.V = np.zeros([sizeX,sizeY]);
        
        
        for i in range(sizeX):
            for j in range(sizeY):
                x = XX[i,j]
                y = YY[i,j]
                dx = x - x0
                dy = y - y0
                r = np.sqrt(dx**2+dy**2)
                self.Vx[i,j] = (lmbda*dx)/(2*np.pi*r**2)
                self.Vy[i,j] = (-lmbda*dy)/(2*np.pi*r**2)
                self.V[i,j]  = np.sqrt(self.Vx[i,j]**2 + self.Vx[i,j]**2)