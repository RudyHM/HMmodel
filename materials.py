# -*- coding: utf-8 -*-


import numpy as np

import library as lib


"""
The file materials includes the definitions of heat and moisture properties
of materials.
    
The following properties are described:
        
'rho': density of the dry material [kg.m-3]  
'k': thermal conductivity [W.m-1.K-1]
'Cs': heat capacity [J.kg-1.K-1]
'delta_l': liquid water permeability [m2]
'delta_v': water vapour permeability [s]
'Xi': hygric capacity [kg.m-3.Pa-1]
"""

# Thermal capacity of the dry material [J.m-3.K-1]
rho_Cp = 1.824e6    
       

def w(T,RH):
    """
    Water content [kg.m-3]   
    Linear interpolation to cover the whole range of relative humidity
    """
    return 146/((1+(8e-8*lib.Pc(T,RH))**(1.6))**0.375)


def delta_l(T,RH):
    """
    Liquid water permeability [s]
    Linear interpolation to cover the whole range of water content
    """            
    return np.exp(-39.2619+0.0704*(w(T,RH)-73)-1.7420e-4*(w(T,RH)-73)**2-2.7953e-6*(w(T,RH)-73)**3-1.1566e-7*(w(T,RH)-73)**4+2.5969e-9*(w(T,RH)-73)**5)
        
  
def delta_v(T,RH):
    """
    Water vapour permeability [s]
    Linear interpolation to cover the whole range of relative humidity
    """
    return lib.Mw/(lib.R*T) * 26.1e-6/200 * (1-w(T,RH)/146)/(0.503*(1-w(T,RH)/146)**2+0.497)

    
def k(T,RH):
    """
    Thermal conductivity [W.m-1.K-1]
    """
    return 1.5 + 15.8/1000*w(T,RH)


    
def Xi(T,RH):
    """
    Hygric capacity (Pc) (dw/dPc) = (dw/dRH*dRH/dPc)
    """     
    return lib.partial_derivative(w,1,[T,RH])*(np.exp(-lib.Pc(T,RH)/(lib.rhoL*lib.Rv*T))/(-lib.rhoL*lib.Rv*T)) 



    
    
    
    
    
    