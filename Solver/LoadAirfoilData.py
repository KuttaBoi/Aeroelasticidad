# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:27:57 2020

@author: qsilv
"""

import numpy as np

def LoadAirfoil(fileName):
    
    hdrlns = 1
    if (fileName == 'nasasc2-0714'):
        hdrlns = 3
    elif (fileName == 's1020'):
        hdrlns = 2
    
    flpth = "C:/Users/qsilv/Documents/Universidad/AAMV/Trabajo Aeroelasticidad/Python/Aeroelasticity/Airfoil_DAT_Selig/"
    flnm = flpth + fileName + ".dat"
    dataBuffer = np.loadtxt(flnm,delimiter=' ', skiprows=hdrlns)
    
    dataX = dataBuffer[:,0]
    dataY = dataBuffer[:,1]
    
    return dataX, dataY