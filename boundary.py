# -*- coding: utf-8 -*-

"""
Boundary conditions
"""

import numpy as np
import xlrd

import library as lib



wb = xlrd.open_workbook('Meteo_Cahors.xlsx')  # opening the Excel file
sh = wb.sheet_by_name(u'Feuil1')                # choosing the sheet   
        
# Reading the file
time = np.asarray(sh.col_values(0))             # time [s]
Text = np.asarray(sh.col_values(1)) + 273.15    # exterior temperature [K]
Tint = np.asarray(sh.col_values(2)) + 273.15    # interior temperature [K]
RHext = np.asarray(sh.col_values(3))            # exterior relative humidity
RHint = np.asarray(sh.col_values(4))            # interior relative humidity
Vs = np.asarray(sh.col_values(5))               # wind speed [m2/s]


# Calculating the capillary pressure 
Pcext = lib.Pc(Text,RHext)  # exterior vapour pressure [Pa]
Pcint = lib.Pc(Tint,RHint)  # interior vapour pressure [Pa]


# Calculating the convective heat transfer coefficient [W/(m2.K)]
hext = 1.7*Vs + 5.1

    
# Creation of a matrice which contains the time, the exterior and interior temperature and vapour pressure
boundary_ini = np.transpose(np.vstack((time, Text, Tint, Pcext, Pcint)))
    
t_tot_ini = int(max(boundary_ini[:,0]))             # total time of the weather data file 
dt_ini = int(boundary_ini[1,0]-boundary_ini[0,0])   # weather data's time step
    
# Implicit method => adding a row in the matrice
boundary_add = np.zeros([1,5])
boundary_add[0,0] = t_tot_ini + dt_ini
for i in range(1,5):
    boundary_add[0,i] = boundary_ini[-1,i]
    
# Creation of a new matrice containing the initial boundary conditions + the adding row
boundary = np.vstack((boundary_ini, boundary_add))   
    
# Parameters regarding the time
dt = int(boundary[1,0]-boundary[0,0])   # time step 
t_tot = int(max(boundary[:,0]))         # total time 


def boundary_interp(dt_sim):
    """
    Linear interpolation of the boundary conditions for the adaptative time step method
    
    Input: simulation time step dt_sim
    """
    Text = np.interp(dt_sim, boundary[:,0], boundary[:,1])
    Tint = np.interp(dt_sim, boundary[:,0], boundary[:,2])
    Pcext = np.interp(dt_sim, boundary[:,0], boundary[:,3])
    Pcint = np.interp(dt_sim, boundary[:,0], boundary[:,4])
    
    return np.vstack((Pcext,Pcint,Text,Tint))
    



    
    
    
    
    
