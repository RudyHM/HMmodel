# -*- coding: utf-8 -*-

"""
Constants and functions
"""

import numpy as np
from scipy.misc import derivative

import mesh


def partial_derivative(func, var=0, point=[]):
    """
    Partial derivative function
    
    Input: function, variable, point
    """
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = 1e-6)
    
    
T_ref = 273.15  # Reference temperature [K] => 0°C
R = 8.314       # Ideal gas constant [J.kg-1.K-1]

# Properties of liquid water
rhoL = 1000.    # Density of liquid water [kg.m-3]
CpL = 4180.     # Heat capacity [J.kg-1.K-1]    
Mw = 18e-3      # Molar mass of water [kg.mol-1]

# Properties of water vapour
# Vapour pressure is defined as a function of both T and RH at the end of this file
CpV = 1850.     # Heat capacity [J.kg-1.K-1] 
Lv = 2.5e6      # Latent heat of evaporation at 0°C [J.kg-1] 
Rv = 462.       # Specific gas constant [J.kg-1.K-1]

# Properties regarding flux - Lewis' analogy to calculate convective mass transfer coefficients
hext = 8.               # Exterior convective heat transfer coefficient [W.m-2.K-1]
hint = 4.               # Interior convective heat transfer coefficient [W.m-2.K-1]
hmext = 7.7e-9 * hext   # Exterior convective mass transfer coefficient [kg.m-2.s-1.Pa-1]
hmint = 7.7e-9 * hint   # Interior convective mass transfer coefficient [kg.m-2.s-1.Pa-1]

    
def Pc(T,RH):
    """
    Capillary pressure [Pa] (choosen positive)  
    Calculated with Gibbs' law 
    
    Input: temperature T [K], relative humidity RH
    """
    return -np.log(RH)*rhoL*Rv*T
   

def Psat(T):
    """
    Water vapour saturation pressure  [Pa]
    Expression from Bristish Standard Institution BS 5250:2002
    
    Input: temperature T [K]
    """       
    return 1000*0.6105*np.exp(17.269*(T-T_ref)/(237.3+(T-T_ref)))
    
    
def dPsat(T):
    """
    Derivative of the water vapour saturation pressure over the temperature
    
    Input: temperature T [K]
    """
    return derivative(Psat, T, dx=1.0)
     
   
def Pv(T,RH):
    """
    Vapour pressure [Pa]
    
    Input: temperature T [K], relative humidity RH
    """
    return RH*Psat(T)

    
def rhoV(T,RH):
    """
    Density of vapour pressure [kg.m-3]
    
    Input: temperature T [K], relative humidity RH
    """
    return Pv(T,RH)/(Rv*T)
    
    
def dPc(T,RH):
    """
    Derivative (first order) of the capillary pressure over x, centered difference
    
    Input: temperature T [K], relative humidity RH
    """
    dPc = np.zeros([mesh.N_tot,1])
     
    for i in range(1,len(T)-1):        
        dPc[i] = (Pc(T,RH)[i+1] - Pc(T,RH)[i-1])/(mesh.dx[i+1]/2+mesh.dx[i]+mesh.dx[i-1]/2)

    return dPc


def d2Pc(T,RH):
    """
    Derivative (second order) of the capillary pressure over x, centered difference
    
    Input: temperature T [K], relative humidity RH
    """    
    d2Pc = np.zeros([mesh.N_tot,1])
    
    for i in range(1,len(T)-1):        
        d2Pc[i] = (Pc(T,RH)[i+1]+Pc(T,RH)[i-1]-2*Pc(T,RH)[i])/(mesh.dx[i+1]/2+mesh.dx[i]+mesh.dx[i-1]/2)**2
    
    return d2Pc


def dT(T):
    """
    Derivative (first order) of the temperature over x, centered difference
    
    Input: temperature T [K]
    """
    dT = np.zeros([mesh.N_tot,1])
     
    for i in range(1,len(T)-1):        
        dT[i] = (T[i+1] - T[i-1])/(mesh.dx[i+1]/2+mesh.dx[i]+mesh.dx[i-1]/2)
         
    return dT
    
    
def d2T(T):
    """
    Derivative (second order) of the temperature over x, centred difference
    
    Input: temperature T [K]
    """    
    d2T = np.zeros([mesh.N_tot,1])
    
    for i in range(1,len(T)-1):        
        d2T[i] = (T[i+1]+T[i-1]-2*T[i])/(mesh.dx[i+1]/2+mesh.dx[i]+mesh.dx[i-1]/2)**2
        
    return d2T



