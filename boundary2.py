# -*- coding: utf-8 -*-

"""
Boundary conditions 
"""

import numpy as np
import xlrd

from library import Pc

div_dt = 60

# Opening the Excel file containing the weather data
wb = xlrd.open_workbook('CL_EN15026.xlsx')
    
# Reading the file
sh = wb.sheet_by_name(u'CL')
time = np.asarray(sh.col_values(0))     # Time [s]
Text = np.asarray(sh.col_values(1))     # Exterior temperature [K]
Tint = np.asarray(sh.col_values(2))     # Interior temperature [K]
RHext = np.asarray(sh.col_values(3))    # Exterior relative humidity
RHint = np.asarray(sh.col_values(4))    # Interior relative humidity


# Vapour pressure is calculated thanks to the temperature and the relative humidity
Pcext = Pc(Text,RHext)                  # Exterior vapour pressure [Pa]
Pcint = Pc(Tint,RHint)                  # Interior vapour pressure [Pa]

# Creation of a matrix which contains the time, the exterior and interior temperature and vapour pressure
boundary_ini = np.transpose(np.vstack((time, Text, Tint, Pcext, Pcint)))

t_tot_ini = int(max(boundary_ini[:,0]))   # Total time of the weather data file 
dt_ini = int(3600)                        # Time step

# Implicit method => adding a row in the matrix 
boundary_plus = np.zeros([1,5])
boundary_plus[0,0] = t_tot_ini + dt_ini
for i in range(1,5):
    boundary_plus[0,i] = boundary_ini[int(t_tot_ini/dt_ini)-1,i]

# Creation of a new matrix containing the initial boundary conditions + the adding row
boundary_sum = np.vstack((boundary_ini, boundary_plus))   

t_tot = int(max(boundary_sum[:,0]))  # New time step


# Simulation's time step
dt = int(dt_ini/div_dt)


def boundary(div_dt):
    """
    Boundary conditions matrix for the simulation     
    
    Input: number which divides the initial time step div_dt
    """
    # Initialization
    boundary = np.zeros([int(div_dt*t_tot/dt_ini),5])
    
    # Copying the values from the matrix boundary_sum (=values from the weather data file)
    for i in range(0,int(t_tot/dt_ini)):
        boundary[i*div_dt] = boundary_sum[i]
                
    # Modification of the first row with the new time step
    for i in range(1,int(div_dt*t_tot/dt_ini)):
        boundary[i,0] = boundary[i-1,0] + dt_ini/div_dt
                
    # Temperature and vapour pressure are interpolated => Only working with constant boundary conditions !!!
    for i in range(1,int(div_dt*t_tot/dt_ini)):
        for j in range(1,5):
            boundary[i,j] = boundary_sum[0,j]
            
    return boundary


boundary = boundary(div_dt)












