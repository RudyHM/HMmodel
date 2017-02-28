# -*- coding: utf-8 -*-

"""
Post-treatment - Plots
"""

import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D  # 3D graphs
from matplotlib import cm, colors  # colors for 3D graphs

import mesh
import boundary as bo


# Results import

StockPc = np.loadtxt('D:/Rudy/These/Modelisation/Validation/EN 15026/Modele/StockPc.txt')
StockT = np.loadtxt('D:/Rudy/These/Modelisation/Validation/EN 15026/Modele/StockT.txt')
StockRH = np.loadtxt('D:/Rudy/These/Modelisation/Validation/EN 15026/Modele/StockRH.txt')
Stockw = np.loadtxt('D:/Rudy/These/Modelisation/Validation/EN 15026/Modele/Stockw.txt')
Stocktsim = np.loadtxt('D:/Rudy/These/Modelisation/Validation/EN 15026/Modele/Stocktsim.txt')


# Creation of storage matrices for 3D plotting
    # The matrices will contain a limited number of values to avoid any memory errors

dt_plot = 300  # Time step for 3D plotting [s]

# Initialization
StockPcplot = np.zeros([mesh.N_tot, int(bo.t_tot/dt_plot)])
StockTplot = np.zeros([mesh.N_tot, int(bo.t_tot/dt_plot)])
StockRHplot = np.zeros([mesh.N_tot, int(bo.t_tot/dt_plot)])
Stockwplot = np.zeros([mesh.N_tot, int(bo.t_tot/dt_plot)])

# First column
StockPcplot[:,0] = StockPc[:,0]
StockTplot[:,0] = StockT[:,0]
StockRHplot[:,0] = StockRH[:,0]
Stockwplot[:,0] = Stockw[:,0]

# Writing
for i in range(1,int(bo.t_tot/dt_plot)):
    StockPcplot[:,i] = StockPc[:,i*int(dt_plot/bo.dt)]
    StockTplot[:,i] = StockT[:,i*int(dt_plot/bo.dt)]
    StockRHplot[:,i] = StockRH[:,i*int(dt_plot/bo.dt)]
    Stockwplot[:,i] = Stockw[:,i*int(dt_plot/bo.dt)]
    

###############################################################################
# 1D plot - line representation at a specific time step

def Pc_line(Mesh_opt):
    """
    Capillary pressure plot at a specific time step
    
    Input: mesh option Mesh_opt
    """
    if Mesh_opt==0:
        Pc_line = plt.plot()
        plt.plot(mesh.Nb_node, StockPc[:,int(bo.t_tot/bo.dt)-1], color='b')
        plt.xlabel('Nodes')
        plt.ylabel('Pc [Pa]')
        
    if Mesh_opt==1 or Mesh_opt==2:
        Pc_line = plt.plot()
        plt.plot(100*mesh.dx_sum, StockPc[:,int(bo.t_tot/bo.dt)-1], color='b')
        plt.xlabel('dx [cm]')
        plt.ylabel('Pc [Pa]')
        
    return Pc_line

    
def T_line(Mesh_opt):
    """
    Temperature plot at a specific time step
    
    Input: mesh option Mesh_opt
    """
    if Mesh_opt==0:
        T_line = plt.plot()
        plt.plot(mesh.Nb_node, StockT[:,int(bo.t_tot/bo.dt)-1], color='r')
        plt.xlabel('Nodes')
        plt.ylabel('T [K]')
        
    if Mesh_opt==1 or Mesh_opt==2:
        T_line = plt.plot()
        plt.plot(100*mesh.dx_sum, StockT[:,int(bo.t_tot/bo.dt)-1], color='r')
        plt.xlabel('dx [cm]')
        plt.ylabel('T [K]')

    return T_line    


def RH_line(Mesh_opt):
    """
    Relative humidity plot at a specific time step
    
    Input: mesh option Mesh_opt
    """
    if Mesh_opt==0:
        RH_line = plt.plot()
        plt.plot(mesh.Nb_node, StockRH[:,int(bo.t_tot/bo.dt)-1], color='g')
        plt.xlabel('Nodes')
        plt.ylabel('RH [-]')

    if Mesh_opt==1 or Mesh_opt==2:
        RH_line = plt.plot()
        plt.plot(100*mesh.dx_sum, StockRH[:,int(bo.t_tot/bo.dt)-1], color='g')
        plt.xlabel('dx [cm]')
        plt.ylabel('RH [-]')
        
    return RH_line


###############################################################################

# 3D plot - Surface representation
    # figure creation
    # base grid creation
    # z axis definition
    # axis labelling
    # graph layout:
        # rstride and cstride: number of lines on the graphs 
        # facecolors: color type
    # figure rotation (elevation, azimuth)
    # plt.gca().invert_xaxis()    # reversing x_axis
    # adding a color bar


def Pc_surf(Mesh_opt):
    """
    Surface representation of the capillary pressure field
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the capillary pressure
    norm = colors.Normalize()    
    
    # Plotting the capillary pressure field
    if Mesh_opt==0:
        Pc_surf = plt.figure()          
        ax = Axes3D(Pc_surf)                        # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = np.arange(0, mesh.N_tot, 1)             # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = StockPcplot[:,:]                        # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing x_axis
        plt.gca().invert_yaxis()    # reversing y_axis
        ax.set_xlabel('time step')
        ax.set_ylabel('Nodes')
        ax.set_zlabel('Pc [Pa]')
        ax.view_init(20, 200)       # figure rotation (elevation, azimuth)
        plt.gca().invert_xaxis()    # reversing x_axis
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('Pc [Pa]', rotation=270)
        
    else:
        Pc_surf = plt.figure()          
        ax = Axes3D(Pc_surf)                        # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = 100*mesh.dx_sum                         # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = StockPcplot[:,:]                        # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing x_axis
        #plt.gca().invert_yaxis()   # reversing y_axis
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('Pc [Pa]')
        ax.view_init(20, -30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('Pc [Pa]', rotation=270)   
        
    return Pc_surf

    
def T_surf(Mesh_opt):
    """
    Surface representation of the temperature field
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the temperature
    norm = colors.Normalize()   
    
    # Plotting the temperature field
    if Mesh_opt==0:
        T_surf = plt.figure()          
        ax = Axes3D(T_surf)                         # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = np.arange(0, mesh.N_tot, 1)             # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = StockTplot[:,:]                         # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        ax.set_xlabel('time step')
        ax.set_ylabel('Nodes')
        ax.set_zlabel('T [K]')
        ax.view_init(20, 200)    # figure rotation (elevation, azimuth)    
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('T [K]', rotation=270)
        
    else:
        T_surf = plt.figure()          
        ax = Axes3D(T_surf)                         # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = 100*mesh.dx_sum                         # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = StockTplot[:,:]                         # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        #plt.gca().invert_yaxis()    # reversing y_axis
        plt.gca().invert_xaxis()    # reversing x_axis
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('T [K]')
        ax.view_init(20, -30)    # figure rotation (elevation, azimuth)    
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('T [K]', rotation=270)    
       
    return T_surf
test = T_surf(1)
    
def RH_surf(Mesh_opt):
    """
    Surface representation of the temperature field
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the relative humidity
    norm = colors.Normalize()  
    
    # Plotting the relative humidity field
    if Mesh_opt==0:
        RH_surf = plt.figure()          
        ax = Axes3D(RH_surf)                        # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = np.arange(0, mesh.N_tot, 1)             # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = StockRHplot[:,:]                        # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        ax.set_xlabel('time step')
        ax.set_ylabel('Nodes')
        ax.set_zlabel('RH [-]')
        ax.view_init(20, 30)    # figure rotation (elevation, azimuth)
            
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('RH [-]', rotation=270)
               
    else:
        RH_surf = plt.figure()          
        ax = Axes3D(RH_surf)                    # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)   # X axis: time representation
        Y = 100*mesh.dx_sum                     # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                 # base grid creation
        Z = StockRHplot[:,:]                    # plot
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_yaxis()    # reversing y_axis  
        #plt.gca().invert_xaxis()    # reversing x_axis    
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('RH [-]')
        ax.view_init(20, 190)    # figure rotation (elevation, azimuth)
        plt.gca().invert_xaxis()    # reversing x_axis
            
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)
        cbar.set_label('RH [-]', rotation=270)
                
    return RH_surf


def w_surf_mat0(Mesh_opt):
    """
    Surface representation of the water content in material 0
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the water content in material 0
    norm = colors.Normalize()  
    
    # Plotting the water content field in material 0
    if Mesh_opt==0:
        w_surf_mat0 = plt.figure()          
        ax = Axes3D(w_surf_mat0)                    # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = np.arange(0, mesh.N_nodes[0], 1)        # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = Stockwplot[0:mesh.N_sum[1],0::]         # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        ax.set_xlabel('time step')
        ax.set_ylabel('Nodes')
        ax.set_zlabel('w_mat0 [kg.m-3]')
        ax.view_init(20, 30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat0 [kg.m-3]', rotation=270)
    
    else:
        w_surf_mat0 = plt.figure()          
        ax = Axes3D(w_surf_mat0)                    # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)  # X axis: time representation
        Y = 100*mesh.dx_sum[0:int(mesh.N_sum[1])]   # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                     # base grid creation
        Z = Stockwplot[0:int(mesh.N_sum[1]),0::]    # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing y_axis       
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('w_mat0 [kg.m-3]')
        ax.view_init(20, -20)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat0 [kg.m-3]', rotation=270)   
       
    return w_surf_mat0


def w_surf_mat1(Mesh_opt):
    """
    Surface representation of the water content in material 1
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the water content in material 1
    norm = colors.Normalize()  
    
    # Plotting the water content field in material 1
    if Mesh_opt==0:
        w_surf_mat1 = plt.figure()          
        ax = Axes3D(w_surf_mat1)                            # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)        # X axis: time representation
        Y = np.arange(0, mesh.N_nodes[1], 1)                # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                             # base grid creation
        Z = Stockwplot[mesh.N_sum[1]:mesh.N_sum[2],0::]     # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        ax.set_xlabel('time step')
        ax.set_ylabel('Nodes')
        ax.set_zlabel('w_mat1 [kg.m-3]')
        ax.view_init(20, 30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat1 [kg.m-3]', rotation=270)

    else:
        w_surf_mat1 = plt.figure()          
        ax = Axes3D(w_surf_mat1)                                    # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)                  # X axis: time representation
        Y = 100*mesh.dx_sum[int(mesh.N_sum[1]):int(mesh.N_sum[2])]  # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                                     # base grid creation
        Z = Stockwplot[int(mesh.N_sum[1]):int(mesh.N_sum[2]),0::]   # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing x_axis   
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('w_mat1 [kg.m-3]')
        ax.view_init(20, -30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat1 [kg.m-3]', rotation=270)

    return w_surf_mat1


def w_surf_mat2(Mesh_opt):
    """
    Surface representation of the water content in material 2
    
    Input: mesh option Mesh_opt
    """
    # Normalizing the colors for the water content in material 2
    norm = colors.Normalize()  
    
    # Plotting the water content field
    if Mesh_opt==0:
        w_surf_mat0 = plt.figure()          
        ax = Axes3D(w_surf_mat0)                            # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot)+1, 1)        # X axis: time representation
        Y = np.arange(0, mesh.N_nodes[2], 1)                # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                             # base grid creation
        Z = Stockwplot[mesh.N_sum[2]:mesh.N_sum[3],0::]     # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing x_axis   
        ax.set_xlabel('time step')
        ax.set_ylabel('Node')
        ax.set_zlabel('w_mat2 [kg.m-3]')
        ax.view_init(20, 30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat2 [kg.m-3]', rotation=270)

    else:
        w_surf_mat0 = plt.figure()          
        ax = Axes3D(w_surf_mat0)                                    # 3D grid initialization
        X = np.arange(0, int(bo.t_tot/dt_plot), 1)                  # X axis: time representation
        Y = 100*mesh.dx_sum[int(mesh.N_sum[2]):int(mesh.N_sum[3])]  # Y axis: space representation
        X, Y = np.meshgrid(X,Y)                                     # base grid creation
        Z = Stockwplot[int(mesh.N_sum[2]):int(mesh.N_sum[3]),0::]   # plot 
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=10, 
                               facecolors=cm.jet(norm(Z)), 
                               linewidth=0,
                               shade=True) 
        plt.gca().invert_xaxis()    # reversing x_axis   
        ax.set_xlabel('time step')
        ax.set_ylabel('dx [cm]')
        ax.set_zlabel('w_mat2 [kg.m-3]')
        ax.view_init(20, 30)    # figure rotation (elevation, azimuth)
        
        # Adding a color bar
        m = cm.ScalarMappable(cmap=surf.cmap, norm=surf.norm)   # Adding a color bar
        m.set_array(Z)
        cbar = plt.colorbar(m)  
        cbar.set_label('w_mat2 [kg.m-3]', rotation=270)
  
    return w_surf_mat2


###############################################################################

def plot_convergence():
    """
    Plotting the convergence of the fields (evolution of the potential over the time at a specific point)
    """
    # Capillary pressure convergence
    Pc_conv = plt.plot()
    plt.plot(StockPc[mesh.N_tot-1,:])
    plt.xlabel('time')
    plt.ylabel('Pcint')
    
    # Temperature convergence
    T_conv = plt.plot()
    plt.plot(StockT[mesh.N_tot-1,:])
    plt.xlabel('time')
    plt.ylabel('Tint')
        
    # Relative humidity convergence
    RH_conv = plt.plot()
    plt.plot(StockRH[mesh.N_tot-1,:])
    plt.xlabel('time')
    plt.ylabel('RHint')
    
    return Pc_conv, T_conv, RH_conv

















