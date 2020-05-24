# -*- coding: utf-8 -*-
"""
******************* UNSTEADY AIRFOIL SOLVER ***************************************
"""

#%% IMPORTAR LIBRERIAS
import numpy as np
import matplotlib.pyplot as plt
import Surface as Surface
from scipy.integrate import solve_ivp
from Grid import Grid
from ElementalSolutions import Vortex
from VORT2D import VORT2D
from unsteady_plotter import PlotAirfoil, PlotWake
from unsteady_model import simulate_airfoil



#%% GEOMETRÍA (Puedes elegir que geometria tiene el perfil)

airfoil = Surface.Parabolic(0,0,1,0,15)
#airfoil = Surface.Airfoil("2032c", 1,0,0)
#airfoil = Surface.Circle(0.1,0,0,16)

#%% PARAMETROS ESTRUCTURALES
b = 1
ah = -0.5
xa = 0.25
ra = 0.5
wh = 1.83
wa = 1
mu = 100
Mass = mu*np.pi*rho*b*b



#%% UNIFORM FLOW
Vinf =01
alfa = 5
alfaR = alfa*(np.pi/180)
wo = 2.1

# DENSIDAD
rho = 1.21

#%% INTEGRACION TEMPORAL
dt = 0.1
t = 0
nTime = 150











































