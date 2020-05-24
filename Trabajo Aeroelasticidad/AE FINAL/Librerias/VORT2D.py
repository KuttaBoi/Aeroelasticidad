import numpy as np

def VORT2D(gamma, x, y, xj, yj):
    rj2 = (x-xj)**2 + (y - yj)**2
    if(rj2<0.00001):
        return 0,0
    u,w = (gamma/(2*np.pi*rj2))*np.array([[0,1],[-1,0]]).dot(np.array([[x-xj],[y-yj]]))
    return u[0],w[0]