import numpy as np
import matplotlib.pyplot as plt

from Panel import Panel

from LoadAirfoilData import LoadAirfoil
from VORT2D import VORT2D

class Surface:
    panels = []    
    c = 1
    trailX = 1.1
    A = 0  
    Ainv = 0
    gamma = []
    dgamma = 0
    
    def __len__(self):
        return len(self.panels)
    
    def PlotSurface(self):
        for panel in self.panels:
            panel.PlotPanel()          
        plt.show()
        
    def PlotNormals(self):
        for panel in self.panels:
            panel.PlotPanelNormals()          
        plt.show()
        
    def PlotPoints(self):
        for panel in self.panels:
            panel.PlotPanelPoints()
        plt.show()
    def DrawShape(self):
        x,y = [p.x0 for p in self.panels],[p.y0 for p in self.panels]
        plt.fill(x,y,'k')

            
        
    def SelfInfluenceMatrix(self):
        numP = len(self.panels)+1
        aij = np.zeros([numP,numP])
        for i, pi in enumerate(self.panels):
            xi = pi.xc
            yi = pi.yc
            for j, pj in enumerate(self.panels):
                xj = pj.xv
                yj = pj.yv
                u,w = VORT2D(1,xi, yi, xj, yj)
                n = pi.normal()
                aij[i,j] = np.dot([u,w],n)
            xj = self.c+0.1
            yj = 0
            u,w = VORT2D(1,xi, yi, xj, yj)
            n = pi.normal()
            aij[i,-1] = np.dot([u,w],n)

        aij[-1,:] = 1
        ainv = np.linalg.inv(aij)
        self.Ainv = ainv
        return aij
    
    def GammaSolve(self,U,V,Wake,h,dh,dtheta,dt):
        b = np.zeros(len(self)+1)
        for i,p in enumerate(self.panels):
            wx, wy = 0, 0
            for v in Wake:
                vx,vy = v.GetSpeedAt(p.xc,p.yc)
                wx = wx + vx 
                wy = wy + vy
                           
            b[i] = -np.array([U + wx - dtheta*h , V + wy + dtheta*(p.xc) - dh]).dot(p.normal())
        b[len(self)] = np.sum([p.gamma for p in self.panels])

        gamma = np.linalg.solve(self.A,b)
        dgamma = (sum(gamma[: -1]) - sum(self.gamma[: -1]))/dt
        self.dgamma = dgamma
        self.gamma = gamma
        for l,p in enumerate(self.panels):
            p.gamma = gamma[l]
        return gamma
    
    def GammaSolve2(self,U,V,Wake):
        b = np.zeros(len(self)+1)
        for i,p in enumerate(self.panels):
            wx = 0
            wy = 0
            for v in Wake:
                vx,vy = v.GetSpeedAt(p.xc,p.yc)
                wx = wx + p.normal().dot(np.array([vx,vy]))   
                           
            b[i] = (U+wx)*np.sin(p.phi)+(V+wy)*np.cos(p.phi)
        b[len(self)] = np.sum([p.gamma for p in self.panels])
        gamma = np.linalg.solve(self.A,b)
        for l,p in enumerate(self.panels):
            p.gamma = gamma[l]
        return gamma
    
    def ForceSolve(self,U,V,dtheta, dh,h,Wake):
        DP = np.zeros(len(self))
        for i, p in enumerate(self.panels): 
            wx = 0
            wy = 0
            for v in Wake:
                vx,vy = v.GetSpeedAt(p.xc,p.yc)
                wx = wx + vx 
                wy = wy + vy
                
            dp = 1.22*(np.array([U + wx - dtheta*h , V + wy + dtheta*(p.xc) - dh]).dot(p.tangential())*p.gamma/p.s + self.dgamma)
            DP[i] = dp
        L = np.sum(DP*[p.s*np.cos(p.phi) for p in self.panels])
        M = -np.sum(DP*[p.s*np.cos(p.phi)*(p.xc-1) for p in self.panels])
        return L,M
    
    def ForceSolve2(self,U,V,Wake):
        DP = np.zeros(len(self))
        for i, p in enumerate(self.panels): 
            wx = 0
            wy = 0
            for v in Wake:
                vx,vy = v.GetSpeedAt(p.xc,p.yc)
                wx = wx + vx 
                wy = wy + vy
                
            dp = 1.22*(U*(p.gamma/p.s) + self.dgamma)
            DP[i] = dp
        L = np.sum(DP)*np.sum([p.s*np.cos(p.phi) for p in self.panels])
        M = -np.sum(DP)*np.sum([p.s*np.cos(p.phi)*(p.xc) for p in self.panels])
        return L,M
        
    

            

        

class Airfoil(Surface):
    size = 1
    
    def __init__(self,filename,size,xOffset,yOffset):
        airfoil = LoadAirfoil(filename)
        XB, YB = airfoil
        XB = XB*size + xOffset
        YB = YB*size + yOffset
        nPan = len(XB) - 1
        panels = [Panel(0,0,0,0) for i in range(len(XB))]      
        
        # La c está normalizada en los archivos y vale 1
        self.c = size
        
        edge = np.zeros(nPan)
        for i in range(nPan):
            edge[i] = (XB[i+1]-XB[i])*(YB[i+1]-YB[i])
        sumEdge = np.sum(edge)
        
        if(sumEdge>0):
            XB = np.flipud(XB)
            YB = np.flipud(YB)
            
        for i in range(nPan):
            panels[i] = Panel(XB[i],YB[i],XB[i+1],YB[i+1])
        self.panels = panels
        self.A = self.SelfInfluenceMatrix()
         

    

class Circle(Surface):
    
    def __init__(self,radius,x0,y0,nIter):
        panels = [Panel(0,0,0,0) for i in range(nIter-1)]
        thetaEnd = (2*np.pi)
        theta = np.linspace(0,thetaEnd,nIter)
        
        XB = np.cos(theta)*radius
        YB = np.sin(theta)*radius
        nPan = len(XB) - 1
        
        edge = np.zeros(nPan)
        for i in range(nPan):
            edge[i] = (XB[i+1]-XB[i])*(YB[i+1]-YB[i])
        sumEdge = np.sum(edge)
        
        if(sumEdge<0):
            XB = np.flipud(XB)
            YB = np.flipud(YB)
            
        for i in range(nPan):
            panels[i] = Panel(XB[i],YB[i],XB[i+1],YB[i+1])
        self.panels = panels
        self.A = self.SelfInfluenceMatrix()
        
        
class Line(Surface):
    
    def __init__(self, x0,y0,x1,y1,nIter):      
        boundX = np.linspace(x0,x1,nIter+1)
        boundY = np.linspace(y0,y1,nIter+1)
        
        self.c = x1-x0
        self.trailX = x1 + 0.2
        
        edge = np.zeros(nIter)
        for i in range(nIter):
            edge[i] = (boundX[i+1]-boundX[i])*(boundY[i+1]-boundY[i])
        sumEdge = np.sum(edge)
        
        if(sumEdge<0):
            boundX = np.flipud(boundX)
            boundY = np.flipud(boundY)
            
        
        for i in range(nIter):
            p = Panel(boundX[i],boundY[i],boundX[i+1],boundY[i+1])
            self.panels.append(p)
        self.A = self.SelfInfluenceMatrix()

class Parabolic(Surface):
    
    def __init__(self, x0,y0,x1,y1,nIter):      
        boundX = np.linspace(x0,x1,nIter+1)
        boundY = np.linspace(y0,y1,nIter+1)
        
        self.c = x1-x0
        self.trailX = x1 + 0.2
        
        boundY = 0.1*4*boundX*(1-boundX)
        
        edge = np.zeros(nIter)
        for i in range(nIter):
            edge[i] = (boundX[i+1]-boundX[i])*(boundY[i+1]-boundY[i])
        sumEdge = np.sum(edge)
        
        if(sumEdge<0):
            boundX = np.flipud(boundX)
            boundY = np.flipud(boundY)
            
        
        for i in range(nIter):
            p = Panel(boundX[i],boundY[i],boundX[i+1],boundY[i+1])
            self.panels.append(p)
        self.A = self.SelfInfluenceMatrix()


        
class SinglePanel(Surface):
    
    def __init__(self, x0,y0,x1,y1):
        panels = [Panel(x0,x1,y0,y1)]
        self.panels = panels
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        