import numpy as np
import matplotlib.pyplot as plt


class Panel:

    # Origen
    x0, y0 = 0, 0
    # Fin
    x1, y1 = 0,0
    # Puntos de control (En el medio)
    xc, yc = 0, 0
    # Puntos de Vortice
    xv, yv = 0, 0
    # Tamaño en x e y
    dx, dy = 0, 0
    s = 1    
    # Ángulos importantes
    phi = 0    
    beta = 0
    # Gamma del Vortex Panel asociado
    gamma = 0
    
    def __init__(self, x0, y0, x1, y1):
        self.x0,self.y0,self.x1,self.y1 = x0,y0,x1,y1
        self.dx, self.dy = x1-x0,y1-y0
        self.s = (self.dx**2+self.dy**2)**0.5
        self.xc,self.yc = self.dx*0.75 + x0, self.dy*0.75 + y0
        self.xv,self.yv = self.dx*0.25 + x0, self.dy*0.25 + y0
        self.phi = -np.arctan2(self.dy,self.dx)

    
    
    def normal(self):
        return np.array([np.sin(self.phi),np.cos(self.phi)])
    def tangential(self):
        return np.array([np.cos(self.phi),-np.sin(self.phi)])
           
    # Ploteamos el panel
    def PlotPanel(self):
        plt.plot([self.x0,self.x1],[self.y0,self.y1],'k', linewidth = 1)
        
    # Ploteamos la normal del panel
    def PlotPanelNormals(self):
        x1 = self.s*self.normal()[0] + self.xc
        y1 = self.s*self.normal()[1] + self.yc 
        plt.plot([self.xc,x1],[self.yc,y1],'r')   
        
    def PlotPanelPoints(self):
        plt.plot(self.x0,self.y0,'ko',markerfacecolor='k',label='Boundary Pts')  
        plt.plot(self.x1,self.y1,'ko',markerfacecolor='k',label='Boundary Pts')  
        
        plt.plot(self.xc,self.yc,'ko',markerfacecolor='b',label='Boundary Pts')  
        plt.plot(self.xv,self.yv,'ko',markerfacecolor='r',label='Boundary Pts')  

        