# -*- coding: utf-8 -*-

"""
Definition of the mesh grid

The multilayered wall is decomposed in N nodes from 0 to N_tot-1
"""

import numpy as np


# Properties regarding the wall
N_L = 1                   # number of layers
etot = np.array([20])     # thickness of each layer [m]                     

Mesh_opt = 1

def Mesh(Mesh_opt):
    
    if Mesh_opt == 0:
        """
        Uniform mesh grid
        """       
        N_nodes = np.array([500])  # number of nodes in each layer
        N_tot = sum(N_nodes[:])    # total number of nodes            
        
        Nb_nodes = np.linspace(0,N_tot-1, num=N_tot, dtype=int)     # numbering of nodes            
        
        # Creating the mesh grid
        dx = np.zeros([N_tot,1])        # initialization

        for i in range(0, N_tot):
            dx[i] = etot/N_tot

        return dx, N_nodes, N_tot, Nb_nodes

        
    if Mesh_opt == 1:    
        """
        Refined mesh grid
        """
        # Refined mesh grid parameters
        u0 = 5e-4   # meshes' thickness at the interface [m] = first term of the geometrical serie
        q = 1.10    # common ration of the geometrical serie
        
        # Initialization of the vector which will contain the number of nodes in each layer     
        N_nodes = int(np.log(1-etot*(1-q)/u0)/np.log(q)-1)            
                   
        # Initializing the vector dx  
        dx = np.zeros([N_nodes,1])      # initialization 
        dx[0] = u0                      # first mesh
            
        for i in range(0, N_nodes-1): 
            dx[i+1] = dx[i]*q

        # Adding a mesh to respect the wall length
        dx_plus = etot-sum(dx[:])
        dx = np.vstack((dx,dx_plus))
        
        # Correcting the number of nodes
        N_nodes = len(dx)

        # Parameters of the mesh grid are reset
        N_tot = N_nodes
        Nb_nodes = np.linspace(0,N_tot-1, num=N_tot)   
        
        return dx, N_nodes, N_tot, Nb_nodes
        
    if Mesh_opt == 2:
        """
        Mesh grid for EN15026
        """
        # Uniform mesh size on the first 0.1m for the water content plot
        dx01 = np.zeros([199,1])
        dx01[:] = 5e-4
        
        # Refined mesh grid parameters
        u0 = 5e-4 # meshes' thickness at the interface [m] = first term of the geometrical serie
        q = 1.10    # common ration of the geometrical serie
        
        # Initialization of the vector which will contain the number of nodes in each layer     
        N_nodes = int(np.log(1-etot*(1-q)/u0)/np.log(q)-1)            
                   
        # Initializing the vector dx  
        dxini = np.zeros([N_nodes,1])      # initialization 
        dxini[0] = u0                      # first mesh
            
        for i in range(0, N_nodes-1): 
            dxini[i+1] = dxini[i]*q
        dx = np.vstack((dx01,dxini))

        # Adding a mesh to respect the wall length
        dx_plus = etot-sum(dxini[:])
        dx = np.vstack((dx,dx_plus))
        
        # Correcting the number of nodes
        N_nodes = len(dx)

        # Parameters of the mesh grid are reset
        N_tot = N_nodes
        Nb_nodes = np.linspace(0,N_tot-1, num=N_tot)   
        
        return dx, N_nodes, N_tot, Nb_nodes
        
# Going out of the function     
dx = Mesh(Mesh_opt)[0]
N_nodes = Mesh(Mesh_opt)[1]
N_tot = Mesh(Mesh_opt)[2]
Nb_nodes = Mesh(Mesh_opt)[3]  


# Defining the cumulative meshes' size vector
dx_sum = np.zeros([N_tot,1])
dx_sum[0] = dx[0]
for i in range(len(dx)):        
    dx_sum[i] = dx_sum[i-1] + dx[i]













