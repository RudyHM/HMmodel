# -*- coding: utf-8 -*-


import numpy as np

import mesh
import materials as mat 
import library as lib


class Matrices(object):
    """
    The class 'Matrices' registers the creation of the matrices for the finite
    difference method. All of the functions below depends on the temperature [K]
    and the relative humidity [-]
    
    'C11': mass conservation coefficient in front of dPc/dt 
    'R11': mass conservation coefficient in front of d²Pc/dx² 
    'R12': mass conservation coefficient in front of d²Pc/dx² 
    'dR11': mass conservation coefficient in front of dPc/dx
    'dR12': mass conservation coefficient in front of dT/dx
    
    'C22': energy conservation coeeficient in front of dT/dt
    'C21': energy conservation coeeficient in front of dPc/dt
    'R22': energy conservation coeeficient in front of d²T/dx²
    'R21': energy conservation coeeficient in front of d²Pc/dx²
    'dR22': energy conservation coeeficient in front of dPc/dx (if k is not a constant)
    'dR21': energy conservation coeeficient in front of dT/dx
    'r1': energy conservation coeeficient in front of (dT/dx)²
    'r2: energy conservation coeeficient in front of dT/dx*dPc/dx
    'R1': linearized energy conservation coefficient in front of dT/dx
    'R2': linearized energy conservation coefficient in front of dPc/dx
    'R3': independant coefficient coming from the linearization of the energy conservation equation
    
    'k0': boundary conditions coefficient
    'k1': boundary conditions coefficient
    
    'Res': resistances matrix (LHS)
    'Capa': capacity matrix (RHS)
    'B': boundary resistances matrix (RHS)
    'K': boundary conditions linearization matrix (RHS)
    'R': linearization of the energy conservation matrix (RHS) 
    """
    
    def __init__(self, T, RH):
        
        # Initialization: temperature and relative humidity are given here
        self.T = T
        self.RH = RH

        
    # Mass conservation coefficients        
    def C11(T,RH):
        """
        Coefficients in front of (dPc/dt)
        
        Input: temperature T [K], relative humidity RH
        """
        return mat.Xi(T,RH)   
        
        
    def R11(T,RH):
        """
        Coefficient in front of (d²Pc/dx²) 
        
        Input: relative humidity RH
        """
        return -(mat.delta_l(T,RH) + mat.delta_v(T,RH) * lib.rhoV(T,RH)/lib.rhoL)

        
    def dR11(R11):
        """
        Coefficient in front of (dPc/dx)
        """
        # Initialization
        dR11 = np.zeros([mesh.N_tot,1])    
        
        # Writing
        for i in range(1,mesh.N_tot-1):
            dR11[i] = (R11[i+1]-R11[i-1])/(mesh.dx[i-1]/2 + mesh.dx[i] + mesh.dx[i+1]/2)
                
        return dR11
        
    def R12(T,RH):
        """
        Coefficient in front of (d²T/dx²)
        """
        return mat.delta_v(T,RH) * (RH*lib.dPsat(T) - lib.Pv(T,RH) * np.log(RH)/T)  

        
    def dR12(R12):
        """
        Coefficient in front of (dT/dx)
        """
        # Initialization
        dR12 = np.zeros([mesh.N_tot,1])    
        
        # Writing
        for i in range(1,mesh.N_tot-1):
            dR12[i] = (R12[i+1]-R12[i-1])/(mesh.dx[i-1]/2 + mesh.dx[i] + mesh.dx[i+1]/2)
            
        return dR12
        
        
    ###########################################################################
    
    # Energy conservation cofficients    
    
    def C22(T,RH):
        """
        Coefficient in front of (dT/dt)
        """
        return mat.rho_Cp + mat.w(T,RH)*lib.CpL
 
       
    def C21(T,RH):
        """
        Coefficient in front of (dPc/dt)
        """
        return ((lib.CpL-lib.CpV)*(T-lib.T_ref) - lib.Lv)*(1./(1-lib.rhoV(T,RH)/lib.rhoL))*mat.Xi(T,RH)
    
    
    def R22(T,RH):
        """
        Coefficient un front of (d²T/dx²)
        """    
        return mat.k(T,RH) 

        
    def dR22(R22):
        """
        Coefficient in front of (dPc/dx)
        
        Input: temperature T [K], relative humidity RH
        """
        # Initialization
        dR22 = np.zeros([mesh.N_tot,1])    
        
        # Writing
        for i in range(1,mesh.N_tot-1):
            dR22[i] = (R22[i+1]-R22[i-1])/(mesh.dx[i-1]/2 + mesh.dx[i] + mesh.dx[i+1]/2)
            
        return dR22      
        
        
    def R21(T,RH):
        """
        Coefficient in front of (d²Pc/dx²)
        """
        return -((lib.CpL-lib.CpV)*(T-lib.T_ref) - lib.Lv) * mat.delta_l(T,RH)

        
    def dR21(R21):
        """
        Coefficient in front of (dPc/dx)
        
        Input: temperature T [K], relative humidity RH
        """
        # Initialization
        dR21 = np.zeros([mesh.N_tot,1])    
        
        # Writing
        for i in range(1,mesh.N_tot-1):
            dR21[i] = (R21[i+1]-R21[i-1])/(mesh.dx[i-1]/2 + mesh.dx[i] + mesh.dx[i+1]/2)
            
        return dR21
              
        
    def r1(T,RH):
        """
        Coefficient in front of (dT/dx)²
        """
        return lib.CpV*mat.delta_v(T,RH)*(RH*lib.dPsat(T)-lib.Pv(T,RH)*np.log(RH)/T)
    
    
    def r2(T,RH):
        """
        Coefficient in front of (dT/dx)*(dPc/dx)
        """
        return -(lib.CpL * mat.delta_l(T,RH) + lib.CpV * mat.delta_v(T,RH) * lib.rhoV(T,RH)/lib.rhoL)
    
    
    def R1(T,RH,r1,r2):
        """
        Coefficient in front of (dT/dx)
        """
        q0 = 2*lib.dT(T)*lib.d2T(T)
        p0 = lib.dPc(T,RH)
        
        return r1*q0 + r2*p0 
    
    
    def R2(T,RH,r2):
        """
        Coefficient in front of (dPc/dx) 
        """
        p1 = lib.dT(T)   
        
        return r2*p1 
    
    
    def R3(T,RH,r1,r2):
        """
        Coefficient coming from the linearization
        """
        q1 = (1 - 2 * lib.d2T(T)) * lib.dT(T)**2
        p2 = - lib.dPc(T,RH) * lib.dT(T) 
        
        return r1*q1 + r2*p2
    
    
    ###########################################################################          
      
    # Matrices       

    def Res(dt, C11, R11, R12, C21, R21, C22, R22, R1, R2):
        """
        Matrix containing the coefficients associated to the driving potential at the instant "j+1"
        """     
        # Initialisation 
        Res1 = np.zeros([mesh.N_tot,mesh.N_tot])  # 1st quarter: Mass conservation
        Res2 = np.zeros([mesh.N_tot,mesh.N_tot])  # 2nd quarter: Influence of heat on water transport
        Res3 = np.zeros([mesh.N_tot,mesh.N_tot])  # 3rd quarter: Influence of water on heat transfer
        Res4 = np.zeros([mesh.N_tot,mesh.N_tot])  # 4th quarter: Energy conservation

        # General case
        for i in range(1,mesh.N_tot-1):
            # Mass conservation
            Res1[i,i+1] = -((R11[i]+R11[i+1])/((mesh.dx[i]/2+mesh.dx[i+1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2))))    # diag+1
            Res1[i,i] = (1./(mesh.dx[i]/2+mesh.dx[i+1]/2+mesh.dx[i]/2+mesh.dx[i-1]/2)*((R11[i+1]+R11[i])/(mesh.dx[i]/2+mesh.dx[i+1]/2) + (R11[i-1]+R11[i])/(mesh.dx[i]/2+mesh.dx[i-1]/2)) + C11[i]/dt)    # diag
            Res1[i,i-1] = -((R11[i]+R11[i-1])/((mesh.dx[i]/2+mesh.dx[i-1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2))))    # diag-1

            # Influence of heat on water transport
            Res2[i,i+1] = -((R12[i]+R12[i+1])/((mesh.dx[i]/2+mesh.dx[i+1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2))))    # diag+1
            Res2[i,i] = (1./(mesh.dx[i]/2+mesh.dx[i+1]/2+mesh.dx[i]/2+mesh.dx[i-1]/2)*((R12[i+1]+R12[i])/(mesh.dx[i]/2+mesh.dx[i+1]/2) + (R12[i-1]+R12[i])/(mesh.dx[i]/2+mesh.dx[i-1]/2)))    # diag
            Res2[i,i-1] = -((R12[i]+R12[i-1])/((mesh.dx[i]/2+mesh.dx[i-1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2))))    # diag-1

            # Influence of water transport on energy conservation
            Res3[i,i+1] = -((R21[i]+R21[i+1])/((mesh.dx[i]/2+mesh.dx[i+1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2)))\
                            + (R2[i+1]+R2[i])/(4*(mesh.dx[i]/2+mesh.dx[i+1]/2)))    # diag+1
            Res3[i,i] = (1./(mesh.dx[i]/2+mesh.dx[i+1]/2+mesh.dx[i]/2+mesh.dx[i-1]/2)*((R21[i+1]+R21[i])/(mesh.dx[i]/2+mesh.dx[i+1]/2) + (R21[i-1]+R21[i])/(mesh.dx[i]/2+mesh.dx[i-1]/2))\
                            + (R2[i+1]+R2[i])/(4*(mesh.dx[i]/2+mesh.dx[i+1]/2)) - (R2[i-1]+R2[i])/(4*(mesh.dx[i]/2+mesh.dx[i-1]/2)) + C21[i]/dt)    # diag
            Res3[i,i-1] = -((R21[i]+R21[i-1])/((mesh.dx[i]/2+mesh.dx[i-1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2)))\
                            - (R2[i-1]+R2[i])/(4*(mesh.dx[i]/2+mesh.dx[i-1]/2)))    # diag-1
            
            # Energy conservation
            Res4[i,i+1] = -((R22[i]+R22[i+1])/((mesh.dx[i]/2+mesh.dx[i+1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2)))\
                            + (R1[i+1]+R1[i])/(4*(mesh.dx[i]/2+mesh.dx[i+1]/2)))    # diag+1
            Res4[i,i] = (1./(mesh.dx[i]/2+mesh.dx[i+1]/2+mesh.dx[i]/2+mesh.dx[i-1]/2)*((R22[i+1]+R22[i])/(mesh.dx[i]/2+mesh.dx[i+1]/2) + (R22[i-1]+R22[i])/(mesh.dx[i]/2+mesh.dx[i-1]/2))\
                        + (R1[i+1]+R1[i])/(4*(mesh.dx[i]/2+mesh.dx[i+1]/2)) - (R1[i-1]+R1[i])/(4*(mesh.dx[i]/2+mesh.dx[i-1]/2)) + C22[i]/dt)    # diag
            Res4[i,i-1] = -((R22[i]+R22[i-1])/((mesh.dx[i]/2+mesh.dx[i-1]/2)*((mesh.dx[i]/2+mesh.dx[i+1]/2)+(mesh.dx[i]/2+mesh.dx[i-1]/2)))\
                        - (R1[i-1]+R1[i])/(4*(mesh.dx[i]/2+mesh.dx[i-1]/2)))    # diag-1

            
        # Boundary conditions

        # Mass conservation
        # Exterior surface
        Res1[0,1] = 0   # diag +1
        Res1[0,0] = 1   # diag
        # Interior surface
        Res1[mesh.N_tot-1,mesh.N_tot-2] = 0     # diag -1
        Res1[mesh.N_tot-1,mesh.N_tot-1] = 1     # diag
  
        # Energy conservation
        # Exterior surface - Neumann boundary conditions
        Res4[0,1] = -2*R22[0]/mesh.dx[0]**2
        Res4[0,0] = 2*R22[0]/mesh.dx[0]**2 + C22[0]/dt + 2*lib.hext/mesh.dx[0]
        # Interior surface
        Res4[mesh.N_tot-1,mesh.N_tot-2] = 0
        Res4[mesh.N_tot-1,mesh.N_tot-1] = 1
        
        # Intermediate concatenation
        Res13 = np.vstack((Res1,Res3))
        Res24 = np.vstack((Res2,Res4))   
         
        return np.hstack((Res13,Res24))
        

    def Capa(dt, C11, C21, C22):    
        """
        Matrix contaning the coefficients associated to the driving potential at the instant "i"
        """
        # Initialisation of each quarter
        Capa1 = np.zeros([mesh.N_tot,mesh.N_tot])     # 1st quarter: Mass conservation
        Capa2 = np.zeros([mesh.N_tot,mesh.N_tot])     # 2nd quarter: Influence of heat on water transport   
        Capa3 = np.zeros([mesh.N_tot,mesh.N_tot])     # 3rd quarter: Influence of water on heat transfer
        Capa4 = np.zeros([mesh.N_tot,mesh.N_tot])     # 4th quarter: Energy conservation
            
        # Writing
        # Mass conservation    
        Capa1 = C11/dt * np.diag(np.ones(mesh.N_tot),0)    
        
        # Influence of water transport on energy conservation   
        Capa3 = C21/dt * np.diag(np.ones(mesh.N_tot),0)    
        
        # Energy conservation
        Capa4 = C22/dt * np.diag(np.ones(mesh.N_tot),0)                
        
        # Dirichlet boundary conditions
        # Mass conservation
        Capa1[0,0] = 0
        Capa1[mesh.N_tot-1,mesh.N_tot-1] = 0
        
        # Influence of water transport on energy conservation
        Capa3[0,0] = 0
        Capa3[mesh.N_tot-1,mesh.N_tot-1] = 0
        
        # Energy conservation
        Capa4[mesh.N_tot-1,mesh.N_tot-1] = 0
        
        # Intermediate concatenation
        Capa13 = np.vstack((Capa1,Capa3))
        Capa24 = np.vstack((Capa2,Capa4))  
        
        return np.hstack((Capa13,Capa24))

        
    def B():
        """
        Matrix associated with the boundary conditions
        """
        # Intialization
        B = np.zeros([2*mesh.N_tot,4])
        
        # Writing
        # Mass conservation boundary condition
        B[0,0] = 1.                 # Exterior boundary
        B[mesh.N_tot-1,1] = 1.      # Interior boundary

        # Energy conservation boundary condition
        B[mesh.N_tot,2] = 2*lib.hext/mesh.dx[0]        # Exterior boundary
        B[2*mesh.N_tot-1,3] = 1.    # Interior boundary

        return B

        
    def R(R3):
        """
        Vector coming from the linearization of the energy conservation equation
        
        Input: temperature T [K], relative humidity RH
        """    
        # Initialization
        R = np.zeros([2*mesh.N_tot,1])
        
        # Writing
        R[mesh.N_tot:2*mesh.N_tot] = -R3
        
        return R
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



