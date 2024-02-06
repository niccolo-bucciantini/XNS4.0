"""
This python routines reads the eos.nb and eos.thermo provided by the COMPOSE library
and relative to a specific EoS and perform a resample of the EoS with uniform spacing in Log10-Log10 scale
both as a function of density, pressure and hentalpy.

BL2 is the EoS by Bombaci & Logotata 2018. (+ one point to make the crust transition a bit smoother)

Plese use with caution and double check with the original COMPOSE table (there might be missing elements!!!!)

If the original COMPOSE table has a coarse sampling the interpolation might not work properly, and the retabulated EoS 
might lead XNS to crash. Make sure the original COMPOSE EoS extends all the way to the low density of the NS Surface.
 
Beware COMPOSE gives in colum 9 the value of : e_tot-1.
"""
import os
import numpy as np
import numpy.ma as ma
import scipy as sp
from scipy.interpolate import interp1d
from scipy.integrate import quad

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob

plt.ion()


# Select the Eos
eosname ='BL2'
# Select the Directory of the EoS
dir=os.getcwd()+'/'+eosname+'/'

# Physical Constants
Msun = 1.99e33
G = 6.6727e-8
c = 3e10
mp = 1.6726219e-24
mn = 1.6749275e-24
me = 9.11e-28
mu = 1.66054e-24
mb = 55.9349375*1.66054e-24/56. + me/2 #from Fe56  
fm = 1.e-13
MeV = 1.6022e-6

# Conversion Factors (1.0005) for SLY and APR (1.0005) for BL2
massfactor = 1.0005

nb2rho = mb/(fm**3)/massfactor  #density to cgs
mevfm32dyne = MeV/(fm**3) #pressure to cgs
eint2ergcm3 = (mn/fm**3)*c**2 #energy to cgs

rhounit = Msun/(G*Msun/c**2)**3
punit = rhounit*c**2

# Read the Data
data1 =[]
data2 =[]
data3 =[]
with open(dir+'eos.nb') as f:
        dummy = int(f.readline())
        npt = int(f.readline())
        for line in f:
            data = line.split()
            data1.append(float(data[0]))

with open(dir+'eos.thermo') as f:
        dummy = f.readline()
        for line in f:
            data = line.split()
            data2.append(float(data[3]))
            data3.append(float(data[9]))

# Set the physical quantities in Geometrized Units
density= np.array(data1)*nb2rho/rhounit            
pressure =  np.array(data2)*np.array(data1)*mevfm32dyne/punit
# Internal & total energy density (beware that COMPOSE gives in column 9 : e_tot-1)
energy = ((1.+np.array(data3))*np.array(data1)*eint2ergcm3-np.array(data1)*nb2rho*c**2)/punit
energyt = ((1.+np.array(data3))*np.array(data1)*eint2ergcm3)/punit

mask_energy =  energy < 0
if mask_energy.sum() > 0:
    print('ATT!!!')
    print('At nb = ',np.array(data1)[np.where(mask_energy)[0][0]],' the internal energy is negative')
    print('Try to change the parameter massfactor, or correct your EoS')
    

# Get the range of the EoS    
minlogrho = np.min(np.log10(density))    
maxlogrho = np.max(np.log10(density))

minlogprs = np.min(np.log10(pressure))    
maxlogprs = np.max(np.log10(pressure))

# Define the resampling points for density in Log10 scale
rhopoints = np.linspace(minlogrho,maxlogrho,1000)
# Define the interpolating function for Log10(P) and Log10(e)
rho2prs = interp1d(np.log10(density),np.log10(pressure))
rho2ein = interp1d(np.log10(density),np.log10(energy))
# Interpolate over the density points
prspoints = rho2prs(rhopoints)
einpoints = rho2ein(rhopoints)

# Compute the local powerlaw index of the EoS (p ~ rho^g; e~ rho^l) in each density interval
prsindex=0.*rhopoints
einindex=0.*rhopoints
prsindex[0:999] = (prspoints[1:1000]-prspoints[0:999])/(rhopoints[1:1000]-rhopoints[0:999])
prsindex[999]=prsindex[998]
einindex[0:999] = (einpoints[1:1000]-einpoints[0:999])/(rhopoints[1:1000]-rhopoints[0:999])
einindex[999]=einindex[998]

# Define dp/drho and dh/drho
dprsdrho = 10**prspoints[0:1000]*prsindex[0:1000]/(10**rhopoints[0:1000])
dhdrho = dprsdrho/(10**rhopoints[0:1000] + 10**prspoints[0:1000] + 10**einpoints[0:1000])

# Define the integral to compute the pseudo-hentalpy
def integrand(x, rhoini, pini, eini, g1, g2):
     # 1/h x dh/drho = 1/(rho +p +e) dp/drho
     return pini*g1*(x/rhoini)**(g1-1.)/rhoini/((x) + pini*(x/rhoini)**g1 + eini*(x/rhoini)**g2)

# Compute the natural log of the hentalpy Ln(h) as a function of density 
I =0*prsindex
LnHent =0*prsindex

pini =  10**prspoints[0]
rhoini = 10**rhopoints[0]
rhofin = 0.
eini = 10**einpoints[0]
g1 = prsindex[0]
g2 = einindex[0]

# Backward integration to set the first point
#print(quad(integrand, rhoini, rhofin, args=(rhoini,pini,eini,g1,g2)))
LnHent[0] = -quad(integrand, rhoini, rhofin, args=(rhoini,pini,eini,g1,g2))[0]
# Integration to set Ln(h)
for i in range(0,999):
    pini =  10**prspoints[i]
    rhoini = 10**rhopoints[i]
    rhofin = 10**rhopoints[i+1]
    eini = 10**einpoints[i]
    g1 = prsindex[i]
    g2 = einindex[i] 

    I[i] = quad(integrand, rhoini, rhofin, args=(rhoini,pini,eini,g1,g2))[0]
    LnHent[i+1]=LnHent[i]+I[i]

# Compute Log10(Ln(h)) and the powelaw index of the EoS (Ln(h) ~ rho^k) in each density interval
entpoints=np.log10(LnHent)
entindex=0.*rhopoints
entindex[0:999] = (entpoints[1:1000]-entpoints[0:999])/(rhopoints[1:1000]-rhopoints[0:999])
entindex[999]=entindex[998]

# Plots in Log10-Log10 scale (check how much in log-log space the trend is linear)
plt.figure('Resampled EoS')
plt.plot(rhopoints,prspoints,label='p')
plt.plot(rhopoints,entpoints,label='Ln[h]')
plt.plot(rhopoints,einpoints,label='e')
plt.xlabel('Log10[density]')
plt.xlabel('Log10[eint], Log10[p], Log10(Ln(h))')
plt.legend()

#---------------------------------------------------------------------------------

# Define the resampling points for pressure in Log10 scale 
prspoints2 = np.linspace(minlogprs,maxlogprs,1000)
# Define the interpolating function For Log10(rho)
prs2rho = interp1d(np.log10(pressure),np.log10(density))
# Interpolate over the pressure points
rhopoints2 = prs2rho(prspoints2)
# Compute the powerlaw index of the EoS (rho ~ p^m) in each pressure interval
rhoindex2=0.*prspoints2
rhoindex2[0:999] = (rhopoints2[1:1000]-rhopoints2[0:999])/(prspoints2[1:1000]-prspoints2[0:999])
rhoindex2[999]=rhoindex2[998]

#---------------------------------------------------------------------------------

# Define the resampling points for Ln(h) in Log10 scale
entpoints3 = np.linspace(entpoints[0],entpoints[999],1000)
# Define the interpolating function For Log10(rho)
ent2rho = interp1d(entpoints,rhopoints)
# Interpolate over the pressure points
rhopoints3 = ent2rho(entpoints3)
# Compute the powerlaw index of the EoS (rho ~ Ln(h)^p) in each pressure interval
rhoindex3=0.*entpoints3
rhoindex3[0:999] = (rhopoints3[1:1000]-rhopoints3[0:999])/(entpoints3[1:1000]-entpoints3[0:999])
rhoindex3[999]=rhoindex3[998]

#---------------------------------------------------------------------------------

# Write the resampled tabulated EoS

npoints=1000
f=open(os.getcwd()+'/'+eosname+"_resampled.dat","w+")
f.write(str(npoints)+'\n')
f.write('\n')
f.write(str('%.15f' %rhopoints[0])+' '+str('%.15f' %rhopoints[999])+'\n')
f.write(str('%.15f' %prspoints[0])+' '+str('%.15f' %prspoints[999])+' '+str('%.15f' %prsindex[0])+' '+str('%.15f' %prsindex[999])+'\n')
f.write(str('%.15f' %einpoints[0])+' '+str('%.15f' %einpoints[999])+' '+str('%.15f' %einindex[0])+' '+str('%.15f' %einindex[999])+'\n')
f.write(str('%.15f' %entpoints[0])+' '+str('%.15f' %entpoints[999])+' '+str('%.15f' %entindex[0])+' '+str('%.15f' %entindex[999])+'\n')
f.write('\n')
# First for uniform Lo10(rho) spacing + powelaw indexes
for i in range(0,1000):
    f.write(str('%.15f' %rhopoints[i])+"\t"+str('%.15f' %prspoints[i])+"\t"+str('%.15f' %einpoints[i])+"\t"+str('%.15f' %entpoints[i])+"\t"+str('%.15f' %prsindex[i])+"\t"+str('%.15f' %einindex[i])+"\t"+str('%.15f' %entindex[i])+'\n')
# First for uniform Lo10(p) spacing + powelaw indexes
f.write('\n')
for i in range(0,1000):
    f.write(str('%.15f' %prspoints2[i])+"\t"+str('%.15f' %rhopoints2[i])+"\t"+str('%.15f' %rhoindex2[i])+'\n')
# First for uniform Lo10(Ln(h)) spacing + powelaw indexes
f.write('\n')
for i in range(0,1000):
    f.write(str('%.15f' %entpoints3[i])+"\t"+str('%.15f' %rhopoints3[i])+"\t"+str('%.15f' %rhoindex3[i])+'\n')
    
f.close()    
