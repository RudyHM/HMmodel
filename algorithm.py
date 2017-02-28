    # -*- coding: utf-8 -*-

"""
Main algorithm

Coupled heat and moisture transfer solved with finite differences method

Validation with the standard EN 15026
"""


###############################################################################
################################### Import ####################################
###############################################################################

# Main libraries
import os 
import numpy as np
import time 

# Work directory    
os.chdir('D:/Users/rbui/Desktop/En cours')    

# Model related imports
import mesh
import boundary as bo
import library as lib
import materials as mat
from method import Matrices as m


###############################################################################
############################### Initialization ################################
###############################################################################

# Definition of the boundary values at the first time step
RHext = bo.RHext[0] 
RHint = bo.RHint[0]
Text = bo.Text[0]
Tint = bo.Tint[0]
Pcext = lib.Pc(Text,RHext)
Pcint = lib.Pc(Tint,RHint)


# Initialization of the initial fields vectors
RHini = np.zeros([mesh.N_tot,1])
Pcini = np.zeros([mesh.N_tot,1])
Tini = np.zeros([mesh.N_tot,1])

initial = 'constant'     # method for the initial conditions

if initial == 'constant':    
    RHini[:] = RHext 
    Pcini[:] = Pcext 
    Tini[:] = Text

if initial == 'linear':
    N = np.linspace(0, mesh.etot, num = mesh.N_tot)   #  abscissa 
    for i in range(0,mesh.N_tot):
        RHini[i] = (RHint-RHext)/mesh.etot * N[i] + RHext
        Tini[i] = (Tint-Text)/mesh.etot * N[i] + Text
        Pcini[i] = (Pcint - Pcext)/mesh.etot * N[i] + Pcext


# Vector [U] which contains Pc (from 0 to N_tot) and T (from N_tot to 2*N_tot-1) => T and Pc are dimensioneless
U = np.zeros([2*mesh.N_tot,1])
U[0:mesh.N_tot,0] = Pcini[:,0]
U[mesh.N_tot:2*mesh.N_tot,0] = Tini[:,0]


# Vectors
Pc = Pcini
T = Tini
RH = RHini


# Time
t_adaptive = 1  # time step option
                    # 0 => dt_sim = dt_ini
                    # 1 => adaptative time step method

# Tolerance for the temporal discretization
Tol = 1e-4


###############################################################################
################################# Resolution ##################################
###############################################################################

# Storage Matrices
StockPc = np.zeros([mesh.N_tot,1])
StockT = np.zeros([mesh.N_tot,1])
StockRH = np.zeros([mesh.N_tot,1])
Stockw = np.zeros([mesh.N_tot,1])

    
# First time step
StockPc = Pc
StockT = T - lib.T_ref
StockRH = RH
Stockw = mat.w(T,RH)


# Resolution
start = time.clock()    # starting the chronometer     
    
if t_adaptive == 0:   
    """
    The simulation time step is equal to the weather data time step
    """
    # Boundary conditions vector
    L = np.zeros([4,1])
    # Implicit method => second value 
    L[0,0] = bo.boundary[1,3]    # Exterior vapour pressure [Pa]
    L[1,0] = bo.boundary[1,4]    # Interior vapour pressure [Pa]
    L[2,0] = bo.boundary[1,1]    # Exterior temperature [K]
    L[3,0] = bo.boundary[1,2]    # Interior temperature [K]
  
    # Loop on every time step
    for k in range(1,int(bo.t_tot/bo.dt)-1):
        
        # Step 1: Creation of the matrices M and N
        M = m.Res(bo.dt, m.C11(T,RH), m.R11(T,RH), m.R12(T,RH), m.C21(T,RH), m.R21(T,RH), m.C22(T,RH), m.R22(T,RH),\
                  m.R1(T,RH,m.r1(T,RH),m.r2(T,RH)), m.R2(T,RH,m.r2(T,RH)), m.dR11(m.R11(T,RH)))
        N = np.dot(m.Capa(bo.dt, m.C11(T,RH), m.C21(T,RH), m.C22(T,RH)),U) + np.dot(m.B(),L) + m.R(m.R3(T,RH,m.r1(T,RH),m.r2(T,RH)))
               
        # Step 2: Solving with LU decomposition
        Uj1 = np.linalg.solve(M,N)
        
        # Step 3: Uptating driving potential vectors
        Pc = Uj1[0:mesh.N_tot]                  # Updating vector Pc
        T = Uj1[mesh.N_tot:2*mesh.N_tot]        # Updating vector T
        RH = np.exp(-Pc/(lib.rhoL*lib.Rv*T))    # Calculating the relative humidity with Gibb's law
        w = mat.w(T,RH)                         # Calculating the water content
        
        # Step 4: Storage of the results 
        StockPc = np.hstack((StockPc,Pc))
        StockT = np.hstack((StockT,T-lib.T_ref))
        StockRH = np.hstack((StockRH,RH))
        Stockw = np.hstack((Stockw,w))
        
        # Step 5: Updating the matrices for the next time step
        # Driving potential vector    
        U = Uj1        
                
        # Boundary conditions matrice
        L[0,0] = bo.boundary[1+k,3]     # Capillary pressure
        L[1,0] = bo.boundary[1+k,4]
        L[2,0] = bo.boundary[1+k,1]     # Temperature
        L[3,0] = bo.boundary[1+k,2]


if t_adaptive == 1:
    """
    The simulation time step is adapted depending on an error calculation and 
    a tolerance factor on the temperature and water content profiles
    """
    # Initialization
    t_sim = 0       # simulation time
    dt_sim = bo.dt  # simulation time step (equal to the weather data time step)
    nb_iter = 0     # number of iteration
    
    Stockdt = np.zeros([1,1])   # storage of the time step
    Stockdt = 0
    StockErr = np.zeros([1,1])  # storage of the error 
    StockErr = 0
  
    while t_sim <= bo.t_tot:     
        
        # Step 1: Coarse calculation from t to t+dt_sim 
        # Boundary conditions vector
        L = bo.boundary_interp(t_sim+dt_sim)
        # Creating matrices M and N
        M = m.Res(dt_sim, m.C11(T,RH), m.R11(T,RH), m.R12(T,RH), m.C21(T,RH), m.R21(T,RH), m.C22(T,RH), m.R22(T,RH),m.R1(T,RH,m.r1(T,RH),m.r2(T,RH)), m.R2(T,RH,m.r2(T,RH)))
        N = np.dot(m.Capa(dt_sim, m.C11(T,RH), m.C21(T,RH), m.C22(T,RH)),U) + np.dot(m.B(),L) + m.R(m.R3(T,RH,m.r1(T,RH),m.r2(T,RH)))
        # Solving with LU decomposition
        U1 = np.linalg.solve(M,N)   # solution of the coarse calculation  
        # Extracting the capillary pressure, temperature and relative humidity
        Pc1 = U1[0:mesh.N_tot]                      # Updating vector Pc
        T1 = U1[mesh.N_tot:2*mesh.N_tot]            # Updating vector T
        RH1 = np.exp(-Pc1/(lib.rhoL*lib.Rv*T1))     # Calculating the relative humidity with Gibb's law
        w1 = mat.w(T1,RH1)                          # Calculating the water content
        

        # Step 2: Refine calculation => the simulation time step is divided by 2
        dt_sim_div = dt_sim/2
        
        # Calculation from t to dt_sim/2
        # Boundary conditions vector
        L = bo.boundary_interp(t_sim + dt_sim_div)
        # Creating matrices M and N
        M = m.Res(dt_sim_div, m.C11(T,RH), m.R11(T,RH), m.R12(T,RH), m.C21(T,RH), m.R21(T,RH), m.C22(T,RH), m.R22(T,RH),\
                  m.R1(T,RH,m.r1(T,RH),m.r2(T,RH)), m.R2(T,RH,m.r2(T,RH)))
        N = np.dot(m.Capa(dt_sim_div, m.C11(T,RH), m.C21(T,RH), m.C22(T,RH)),U) + np.dot(m.B(),L) + m.R(m.R3(T,RH,m.r1(T,RH),m.r2(T,RH)))
        # Solving with LU decompositon
        Uinterm = np.linalg.solve(M,N)  # intermediate solution 
        # Extracting the capillary pressure, temperature and relative humidity
        Pcinterm = Uinterm[0:mesh.N_tot]                        # Updating vector Pc
        Tinterm = Uinterm[mesh.N_tot:2*mesh.N_tot]              # Updating vector T
        RHinterm = np.exp(-Pcinterm/(lib.rhoL*lib.Rv*Tinterm))  # Calculating the relative humidity with Gibb's law   
        
        # Calculation from dt_sim/2 to t+dt_sim using Uinterm 
        L = bo.boundary_interp(t_sim + 2*dt_sim_div)    # Updating the boundary conditions vector
        # Creating matrices M and N
        M = m.Res(dt_sim_div, m.C11(Tinterm,RHinterm), m.R11(Tinterm,RHinterm), m.R12(Tinterm,RHinterm), m.C21(Tinterm,RHinterm), m.R21(Tinterm,RHinterm), m.C22(Tinterm,RHinterm),\
                  m.R22(Tinterm,RHinterm), m.R1(Tinterm,RHinterm,m.r1(Tinterm,RHinterm),m.r2(Tinterm,RHinterm)), m.R2(Tinterm,RHinterm,m.r2(Tinterm,RH)))
        N = np.dot(m.Capa(dt_sim_div, m.C11(Tinterm,RHinterm), m.C21(Tinterm,RHinterm), m.C22(Tinterm,RHinterm)),Uinterm) + np.dot(m.B(),L)\
            + m.R(m.R3(Tinterm,RHinterm,m.r1(Tinterm,RHinterm),m.r2(Tinterm,RHinterm)))
        # Solving with LU decomposition
        U2 = np.linalg.solve(M,N)   # solution of the refine calculation  
        # Extracting the capillary pressure, temperature and relative humidity
        Pc = U2[0:mesh.N_tot]                   # Updating vector Pc
        T = U2[mesh.N_tot:2*mesh.N_tot]         # Updating vector T
        RH = np.exp(-Pc/(lib.rhoL*lib.Rv*T))    # Calculating the relative humidity with Gibb's law
        w = mat.w(T,RH)                         # Calculating the water content
      

        # Step 3: Error calculation between U1 and U2
        # Errors on the temperature and capillary pressure fields are calculated independently
        # The maximum error is selected
        ErrPc = max(abs((U1[0:mesh.N_tot]-U2[0:mesh.N_tot])/U2[0:mesh.N_tot]))
        ErrT = max(abs((U1[mesh.N_tot:2*mesh.N_tot]-U2[mesh.N_tot:2*mesh.N_tot])/U2[mesh.N_tot:2*mesh.N_tot]))
        Err = max(ErrPc,ErrT)

        # Step 5: Updating the time step
        if Err > Tol:
            """
            The error is too big
            """
            nb_iter = nb_iter+1     # keeping track of the number of iteration
            dt_sim = dt_sim_div
        
        else:
            """
            The tolerance criterion is respected => the optimal time step is used
            """     
            # Updating the driving potential vector
            U = U2  
            
            # Storing the results
            StockPc = np.hstack((StockPc,Pc))
            StockT = np.hstack((StockT,T-lib.T_ref))
            StockRH = np.hstack((StockRH,RH))
            Stockw = np.hstack((Stockw,w))        
            Stockdt = np.hstack((Stockdt,dt_sim))
            StockErr = np.hstack((StockErr,Err))
            
            nb_iter = nb_iter+1     # keeping track of the number of iteration
            t_sim = t_sim + dt_sim  # updating the simulation time

            dt_opt = 0.9*dt_sim*(Tol/Err)**(1./2)   # ponderated by 0.9 for a safety
            """
            if dt_opt > 43200:  # upper limit: 12h (43200s)
                dt_opt = 43200
            """
            if t_sim > bo.t_tot:
                """
                Possibility that the total simulation time exceed the total time 
                of the data file. If it happens, the linear interpolation is 
                still done by considering constant boundary conditions and equal
                to the last step
                """          
                dt_sim = dt_sim - (t_sim - bo.t_tot)
                Stockdt[-1] = dt_sim    # correction of the last time step 

            else:
                dt_sim = dt_opt

    t_sim_tot = sum(Stockdt[:])     # checking if the simulation time match the total time of the boundary conditions file
                
    # Cumulative time step matrice
    Stocktsim = np.zeros([len(Stockdt),1])  # in seconds
    Stocktsim[0] = Stockdt[0]
    for i in range(0,len(Stockdt)-1):
        Stocktsim[i+1] = Stocktsim[i] + Stockdt[i+1]

    Stocktsimd = Stocktsim/(24*3600)    # conversion in days
                    
elapsed = (time.clock() - start)    # end of the chronometer

"""
# Saving data in text files
np.savetxt('StockPc.txt', StockPc)
np.savetxt('StockT.txt', StockT)
np.savetxt('StockRH.txt', StockRH)
np.savetxt('Stockw.txt', Stockw)
np.savetxt('StockRH.txt', StockRH)
if t_adaptive == 1:
    np.savetxt('Stocktsim.txt', Stocktsim)
"""















