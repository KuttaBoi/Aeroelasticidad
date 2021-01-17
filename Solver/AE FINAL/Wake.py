# -*- coding: utf-8 -*-
"""
Created on Sun May 24 10:23:30 2020

@author: qsilv
"""
import numpy as np
from VORT2D import VORT2D

#%% ESTELA
class Wake:
    Wake = []
    RemoveSmallVortices = 1
    airfoil = 0
    
    def __init__(self,airfoil):
        self.airfoil = airfoil
    
    def UpdateWakePosition(self, Ut,Wt,dt):
        for w in Wake:
            #Eliminamos vortices chiquitos?
            if(self.RemoveSmallVortices and (np.abs(w.gamma)<0.001 or w.x > 20)):
                Wake.remove(w)
                break
            vp, wp = 0,0 
            for p in self.airfoil.panels:
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