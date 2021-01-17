# ------------ GRID ----------------

import numpy as np
import matplotlib.pyplot as plt


class Grid:
    
    XX = []
    YY = []
    X = []
    Y = []
    sizeX = 0
    sizeY = 0
    numX = 0
    numY = 0
    
    def __init__(self, sizeX, sizeY, numX, numY):
        
        self.X = np.linspace(sizeX[0],sizeX[1],numX)
        self.Y = np.linspace(sizeY[0],sizeY[1],numY)
        self.XX, self.YY = np.meshgrid(self.X,self.Y)   
        self.sizeX = sizeX
        self.sizeY = sizeY
        self.numX = numX
        self.numY = numY
        
    def PlotStream(self,Fx,Fy):
        numbOfStreamLines = 50
        slX = np.ones(numbOfStreamLines)
        slY = np.linspace(self.sizeY[0],self.sizeY[1],numbOfStreamLines)
        slXY = np.vstack((slX.T,slY.T)).T
        plt.streamplot(self.XX,self.YY,Fx,Fy,density = 10, linewidth=0.4,color='g',arrowstyle='-',start_points=slXY)
        