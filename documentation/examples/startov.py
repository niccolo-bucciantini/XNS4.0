#!/usr/bin/env python
# coding: utf-8

# Profiles of the main quantities of the XNS 4.0 code</center>
#
# In this notebook we plot the radial profile of the main outputs of the XNSmod code.

# In[9]:

import warnings
#warnings.simplefilter("error")
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from decimal import Decimal
from scipy.interpolate import splprep, splev

# Conversion from code to physical units in cgs
rhounit = 6.2031e17 
bunit   = 8.3758e19
runit   = 1.47459

# Range of the Plot in km (for other settings see below)
RangePlot =18
# Plot variable -  Set True for saving to png file
PLOT = False

# =============================================================================
# CHECK
# =============================================================================

cwd       = os.getcwd()
gridpath  = os.path.join(cwd,'Grid.dat')
datapath = os.path.join(cwd,'TOVINIMOD_PROFILES.dat')

# Check for correct working directory
verify = np.logical_and(os.path.exists(datapath),os.path.exists(gridpath))
if np.logical_not(verify):
    print('ATT! please run this file in the directory containing Grid.dat & TOVINIMOD_PROFILES.dat')
    sys.exit()
 

# =============================================================================
# DICTIONARIES
# =============================================================================

#Colorbar orientation
orientation = {'v': 'vertical', 'h': 'horizontal'}

# =============================================================================
# SETTING - TO BE CUSTROMIZED BY USER
# =============================================================================

#Figure size
xsize = 8
ysize = 8

# Plot Range (same in X and y)
Xmin = 0.0
Xmax =  RangePlot

# Font size
FontSizeAxis  = 10
FontSizeTitle = 10
lbs = 15
fs  = lbs

# "Tick frequency" on x and y axis
tfx = 5 #x-axis
tfy = 5 #y-axis

# Plot variable -  Set True for saving to png filke
PLOT = False

# Contour lines color and linestyles in contour plot
rgb_r = 16  #rgb red
rgb_g = 52  #rgb green
rgb_b = 166 #rgb blue
COLORS = (rgb_r/255,rgb_g/255,rgb_b/255)
LINESTYLES = 'solid'

# Color map and transparency in contour plot
CMAP = 'jet'
ALPHAC = 1.0 #(0<= alpha <= 1)
fr = 0.043
pd = -0.1

# Settings for streamplot
rgb_r_sp = 16 #rgb red
rgb_g_sp = 52 #rgb green
rgb_b_sp = 166 #rgb blue
COLOR_SP = (rgb_r_sp/255,rgb_g_sp/255,rgb_b_sp/255)
ARROWSIZE = 1.5
LINEWIDTH = 1.
DENSITY = 150


# ### Grid, Main and Read functions

# =============================================================================
# GRID FUNCTION
# Read Grid.dat and set the grid for plotting
# =============================================================================
def grid():
    f = np.fromfile(os.path.join(cwd,'Grid.dat'),count=7,sep=' ')
    
    NTH = int(f[0])
    NR = int(f[1])

    data = np.genfromtxt(os.path.join(cwd,'Grid.dat'), skip_header=1, usecols=0, unpack=True)
    
    TH = data[0:NTH]
    R  = data[NTH:NTH+NR]
    Rmin = R.min()
    Rmax = R.max()
    rad = np.insert(R,0,-R[0])
    rad = np.insert(rad,NR+1,R[NR-1]+0.5*(R[NR-1]-R[NR-2]))
    radius = 0.5*(rad[1:NR+2]+rad[0:NR+1])
    theta  = np.linspace(0,np.pi,NTH+1)

    r2,th2 = np.meshgrid(radius,theta)
    XC = np.outer(np.sin(theta),radius)
    YC = np.outer(np.cos(theta),radius)
    return NTH,NR,R,TH,XC,YC,Rmin,Rmax
  
# =============================================================================
# MAIN FUNCTION
# 
# =============================================================================

def main():
    #Grid
    NTH, NR, R, TH, XC, YC, Rmin, Rmax = grid()

    IZ  = np.genfromtxt(datapath, skip_header=1, usecols=0, unpack=True)
    RHO = np.genfromtxt(datapath, skip_header=1, usecols=1, unpack=True)
    PRS = np.genfromtxt(datapath, skip_header=1, usecols=2, unpack=True)
    NU  = np.genfromtxt(datapath, skip_header=1, usecols=3, unpack=True)
    MU  = np.genfromtxt(datapath, skip_header=1, usecols=4, unpack=True)
    CHI = np.genfromtxt(datapath, skip_header=1, usecols=5, unpack=True)

    #Figure size
    fig, ax = plt.subplots(1,1,figsize=(xsize,ysize))
    plt.xlabel('R $[$km$]$', fontsize=lbs)
    plt.subplots_adjust(bottom=0.15)
    ax = fig.gca()
    ax.set_xlim(Xmin,Xmax)

    Field = input('Density  [r], Pressure [p], nu [nu], mu [mu], chi [ch]: ')
  
                
    if Field == 'r':
        print('RHO max: ', '{:.2e}'.format(Decimal(RHO.max()*rhounit)),'g/cm^3')
        filename='Density_tov.png'
        #Colored map of RHO
        # PLTF=np.log10(RHO+1.E-10)
        PLTF=RHO*rhounit/1.0e14
        title = '$\\rho$ $[10^{14}$g/cm$^3$$]$'
        plt.ylabel('$\\rho$ $[10^{14}$g/cm$^3$$]$', fontsize=lbs)
        ax.set_ylim(0,1.1*PLTF.max())

    if Field == 'p':
        print('PRS max: ', '{:.2e}'.format(Decimal(PRS.max()*rhounit)),'g/cm^3')
        filename='Pressure_tov.png'
        #Colored map of PRS
        # PLTF=np.log10(PRS+1.E-14)
        PLTF=PRS*rhounit/1.0e14
        title = '$prs/c^2$ $[10^{14}$g/cm$^3$$]$'
        plt.ylabel('$prs/c^2$ $[10^{14}$g/cm$^3$$]$', fontsize=lbs)
        ax.set_ylim(0,1.1*PLTF.max())

    
    if Field == 'nu':
        print('VPHI max: ', '{:.2e}'.format(Decimal(VPHI.max())),'[c]')
        filename='Nu_tov.png'
        #Colored map of NU
        PLTF= NU
        title = '$\\nu$'
        plt.ylabel(title, fontsize=lbs)
        ax.set_ylim(1.1*PLTF.min(),1.1*PLTF.max())

    if Field == 'mu':
        print('Alpha min: ', '{:.2e}'.format(Decimal(ALPHA.min())),)
        filename='Mu_tov.png'
        #Colored map of MU
        PLTF= MU
        title = '$\\mu$'
        plt.ylabel(title, fontsize=lbs)
        ax.set_ylim(1.1*PLTF.min(),1.1*PLTF.max())

    if Field == 'ch':
        print('Psi max: ', '{:.2e}'.format(Decimal(PSI.max())),)
        filename='Chi_tov.png'
        #Colored map of CHI
        PLTF= CHI
        title = '$\\chi$'
        plt.ylabel(title, fontsize=lbs)
        ax.set_ylim(0,1.1*PLTF.max())


    plt.plot(R[:NR-1]*runit,PLTF[:],color='blue',linewidth=2)

   #Showing
    if (PLOT):
        fig.savefig(os.path.join(cwd,filename), format='png', dpi=70, bbox_inches='tight')

    plt.show()

    
# =============================================================================
# RUN FILE AS MAIN
# =============================================================================
if __name__ == "__main__":
    main()
