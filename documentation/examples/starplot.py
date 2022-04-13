# This notebook is use for 2D plot the main outputs of the XNS code

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
RangePlot = 30
# Plot variable -  Set True for saving to png file
PLOT = True

# =============================================================================
# CHECK
# =============================================================================

cwd       = os.getcwd()
gridpath  = os.path.join(cwd,'Grid.dat')
hydropath = os.path.join(cwd,'Hydroeq.dat')
hydromagpath = os.path.join(cwd,'Hydroeq_mag.dat')
surfpath  = os.path.join(cwd,'Surf.dat')

# Check for correct working directory
verify = np.logical_and(os.path.exists(hydropath),os.path.exists(gridpath))
if np.logical_not(verify):
    print('ATT! please run this file in the directory containing Grid.dat & Hydroeq.dat')
    sys.exit()
verify = os.path.exists(surfpath)
if np.logical_not(verify):
    print('ATT! please run this file in the directory containing Surf.dat')
    sys.exit()    
magnetized = os.path.exists(hydromagpath) 
if np.logical_not(magnetized):
    print('The directory does not contain Hydroeq_mag.dat')
    print('This is not a Magnetized Model')
    print('Plot is only for Hydro amnd Metric Variables')
 

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
Xmin = -RangePlot
Xmax =  RangePlot
Ymin = -RangePlot
Ymax =  RangePlot

# Label x-axis (or both axes if contour plot)
XAxisStr = 'x-axis'
YAxisStr = 'z-axis'

# Font size
FontSizeAxis  = 10
FontSizeTitle = 10
lbs = 15
fs  = lbs

# "Tick frequency" on x and y axis
tfx = 5 #x-axis
tfy = 5 #y-axis

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
# COVTERM FUNCTION
# Compute the 3-metric terms gamma_rr, gamma_tt, gamma_pp
# =============================================================================
def covterm(psi,R,TH,ascal2):
    GCOVR = psi**4*ascal2
    GCOVT = psi**4*np.square(R)*ascal2
    GCOVP = psi**4*np.square(R)*np.square(np.sin(TH))*ascal2
    return GCOVR,GCOVT,GCOVP

# =============================================================================
# READ FUNCTION
# Read Hydroeq.dat, Hydroeq_mag.dat and Surf.dat
# =============================================================================
def read(NTH,NR,R,TH,LINES=False):
    # density rho
    RHO = np.genfromtxt(hydropath, skip_header=1, usecols=0, unpack=True).reshape(NTH,NR)
    # pressure prs
    PRS = np.genfromtxt(hydropath, skip_header=1, usecols=1, unpack=True).reshape(NTH,NR)
    # psi function
    PSI = np.genfromtxt(hydropath, skip_header=1, usecols=2, unpack=True).reshape(NTH,NR)
    # v^phi function
    v3 = np.genfromtxt(hydropath, skip_header=1, usecols=3, unpack=True).reshape(NTH,NR)
    # alpha
    ALPHA = np.genfromtxt(hydropath, skip_header=1, usecols=4, unpack=True).reshape(NTH,NR)
    # beta^phi
    BETA  = np.genfromtxt(hydropath, skip_header=1, usecols=5, unpack=True).reshape(NTH,NR)
    # chi function
    CHI = np.genfromtxt(hydropath, skip_header=1, usecols=6, unpack=True).reshape(NTH,NR)
    ascal2=np.exp(-2.0E-4*CHI*2-6*CHI**2)
    # dchi/dr function
    QSCALR = np.genfromtxt(hydropath, skip_header=1, usecols=7, unpack=True).reshape(NTH,NR)
    # dchi/dtheta function
    QSCALT = np.genfromtxt(hydropath, skip_header=1, usecols=8, unpack=True).reshape(NTH,NR)
    
    # stellar surface
    SURF = np.genfromtxt(surfpath, usecols=0, unpack=True)

    
    VPHI  = np.zeros((NTH,NR))
    B2    = np.zeros((NTH,NR))
    B2POL = np.zeros((NTH,NR))
    B2TOR = np.zeros((NTH,NR))
    E2POL = np.zeros((NTH,NR))
    J2POL = np.zeros((NTH,NR))
    J2TOR = np.zeros((NTH,NR))

    if magnetized:
        #Magnetic Field (Hydroeq_mag file)
        dataMagB = np.genfromtxt(hydromagpath, skip_header=1, usecols=(0,1,2), unpack=True)
        APHI = np.genfromtxt(hydropath, skip_header=1, usecols=3, unpack=True).reshape(NTH,NR)
        dataMagE = np.genfromtxt(hydromagpath, skip_header=1, usecols=(4,5,6), unpack=True)
        ATIM = np.genfromtxt(hydropath, skip_header=1, usecols=7, unpack=True).reshape(NTH,NR)
        dataMagJ = np.genfromtxt(hydromagpath, skip_header=1, usecols=(8,9,10), unpack=True)
        
        b3    = dataMagB[0].reshape(NTH,NR)
        bpolr = dataMagB[1].reshape(NTH,NR)
        bpolt = dataMagB[2].reshape(NTH,NR)

        e3    = dataMagE[0].reshape(NTH,NR)
        epolr = dataMagE[1].reshape(NTH,NR)
        epolt = dataMagE[2].reshape(NTH,NR)

        j3    = dataMagJ[0].reshape(NTH,NR)
        jpolr = dataMagJ[1].reshape(NTH,NR)
        jpolt = dataMagJ[2].reshape(NTH,NR)
    
    for i in range(0,NTH):
        for k in range(0,NR):
            GCOVR,GCOVT,GCOVP = covterm(PSI[i,k],R[k],TH[i],ascal2[i,k])
            VPHI[i,k] = np.square(v3[i,k])*GCOVP
            if(magnetized):
                B2POL[i,k] = np.square(bpolr[i,k])*GCOVR + np.square(bpolt[i,k])*GCOVT
                B2TOR[i,k] = np.square(b3[i,k])*GCOVP
                B2[i,k]    = B2POL[i,k] + B2TOR[i,k]

                E2POL[i,k] = np.square(epolr[i,k])/GCOVR + np.square(epolt[i,k])/GCOVT

                J2POL[i,k] = np.square(jpolr[i,k])*GCOVR + np.square(jpolt[i,k])*GCOVT
                J2TOR[i,k] = np.square(j3[i,k])*GCOVP

                
        
    BTOT = np.sqrt(B2)   
    BPOL = np.sqrt(B2POL)
    BTOR = np.sqrt(B2TOR)

    EPOL = np.sqrt(E2POL) 

    JPOL = np.sqrt(J2POL)
    JTOR = np.sqrt(J2TOR)
    
    if LINES == False:
        return RHO, PRS, PSI, VPHI, ALPHA, BETA, CHI, BTOT, BTOR, BPOL, EPOL, JPOL, JTOR, SURF
    elif LINES == True:
        return bpolr,bpolt,PSI,BTOT,CHI,SURF,RHO

# =============================================================================
# FIELDLINES FUNCTION
# Draw the filedlines of a poloidal field
# =============================================================================    

def FieldLines(NTH,NR,R,TH,Rmin,Rmax):   
    bpolr,bpolt,psi,B,chi,surf,rho = read(NTH,NR,R,TH,LINES=True)

    Rint=int(Rmax/2)
    brr = np.zeros((Rint,Rint))
    btt = np.zeros((Rint,Rint))
    
    x = np.linspace(-Rint/np.sqrt(2),Rint/np.sqrt(2),Rint)
    z = np.linspace(-Rint/np.sqrt(2),Rint/np.sqrt(2),Rint)
    
    for i in range(0,Rint):
        for k in range(0,Rint):
            xpt = x[i]
            zpt = z[k]
            rpt = np.sqrt(xpt**2+zpt**2)
            thpt = np.arccos(zpt/(rpt+1.e-12))
            rint = np.interp(rpt,R,np.linspace(1,NR,NR))
            rint = rint.astype(int)
            thint  = (thpt/np.pi*NTH).astype(int)
            alphas = np.exp(-0.0002*chi[thint,rint]-6/2.*chi[thint,rint]**2)
            brloc  = bpolr[thint,rint]
            bthloc = bpolt[thint,rint]
            psiloc = psi[thint,rint]                    
            bxx = np.sign(xpt)*(brloc*np.sin(thpt) + bthloc*rpt*np.cos(thpt))
            bzz = brloc*np.cos(thpt) - bthloc*rpt*np.sin(thpt)
            brr[k,i] = bxx*psiloc**2*alphas
            btt[k,i] = bzz*psiloc**2*alphas
 
    en = np.sqrt(brr**2 + btt**2)
    lw = 2*en / en.max() +.5

    xstart = np.linspace(0.5,10,15)
    ystart = np.linspace(-0,0,15)
    startp = np.vstack([xstart,ystart]).T
    
    strm = plt.streamplot(x*runit,z*runit,brr*bunit,btt*bunit,color='white',linewidth=0.5,arrowsize=0.00001,density=30,start_points=startp)
    
    xstart = np.linspace(-10,-0.5,15)
    ystart = np.linspace(-0,0,15)
    startp = np.vstack([xstart,ystart]).T
    
    strm = plt.streamplot(x*runit,z*runit,brr*bunit,btt*bunit,color='white',linewidth=0.5,arrowsize=0.00001,density=30,start_points=startp)
    
# =============================================================================
# MAIN FUNCTION
# 
# =============================================================================

def main():
    #Grid
    NTH, NR, R, TH, XC, YC, Rmin, Rmax = grid()

    #Magnetic and Scalar Fields
    RHO, PRS, PSI, VPHI, ALPHA, BETA, CHI, BTOT, BTOR, BPOL, EPOL, JPOL, JTOR, SURF = read(NTH,NR,R,TH)

    BpolMaxRef = np.max(np.abs(BPOL))
    VelMaxRef = np.max(np.abs(VPHI))
    #Figure size
    fig, ax = plt.subplots(1,1,figsize=(xsize,ysize))
    plt.ylabel('z $[$km$]$', fontsize=lbs)
    plt.xlabel('x $[$km$]$', fontsize=lbs)
    plt.subplots_adjust(bottom=0.15)
    ax = fig.gca()
    ax.set_aspect('equal')
    ax.set_xlim(Xmin,Xmax)
    ax.set_ylim(Ymin,Ymax)


    print('')
    if magnetized:
        VarSet = input('Fluid Vars  [f], Metric Vars [g], or Magnetic Vars [m]: ')
    else:
        VarSet = input('Fluid Vars  [f], Metric Vars [g]: ')
        
    if VarSet == 'f':
        Field = input('Density [r], Pressure [p], Velocity [v]: ')
    if VarSet == 'g':
        Field = input('alpha [al], psi [ps], beta [bs], chi[ch]: ')
    if magnetized:
        if VarSet == 'm':
            if np.logical_and(BpolMaxRef >= 1.e-12,  VelMaxRef<=1.e-12):
                print('This non rotating model either with poloidal or twisted-torus magnetic field')
                print('No electric field')
                Field = input('B [bb], Bpol [bp], Btor [bt], Jpol[jp], Jtor[jt]: ')
            elif np.logical_and(BpolMaxRef >= 1.e-12,  VelMaxRef>=1.e-12):
                print('This rotating model with purely poloidal magnetic field')
                print('Currents are not provided in this case')
                Field = input('B [bb], Bpol [bp], Epol[ep]: ')
            else:
                print('This model has a purely toroidal magnetic field')
                print('Currents are not provided in this case')
                Field = input('B [bb], Btor [bt]: ')
                
    if Field == 'r':
        print('RHO max: ', '{:.2e}'.format(Decimal(RHO.max()*rhounit)),'g/cm^3')
        filename='Density.png'
        #Colored map of RHO
        # PLTF=np.log10(RHO+1.E-10)
        PLTF=RHO*rhounit/1.0e14
        title = '$\\rho$ $[10^{14}$g/cm$^3$$]$' 

    if Field == 'p':
        print('PRS max: ', '{:.2e}'.format(Decimal(PRS.max()*rhounit)),'g/cm^3')
        filename='Pressure.png'
        #Colored map of PRS
        # PLTF=np.log10(PRS+1.E-14)
        PLTF=PRS*rhounit/1.0e14
        title = '$prs/c^2$ $[10^{14}$g/cm$^3$$]$'
    
    if Field == 'v':
        print('VPHI max: ', '{:.2e}'.format(Decimal(VPHI.max())),'[c]')
        filename='Velocity.png'
        #Colored map of VPHI
        PLTF= VPHI
        title = '$v/c$'

    if Field == 'al':
        print('Alpha min: ', '{:.2e}'.format(Decimal(ALPHA.min())),)
        filename='Lapse.png'
        #Colored map of ALPHA
        PLTF= ALPHA
        title = '$\\alpha$'

    if Field == 'ps':
        print('Psi max: ', '{:.2e}'.format(Decimal(PSI.max())),)
        filename='Conformal.png'
        #Colored map of PSI
        PLTF= PSI
        title = '$\\psi$'

    if Field == 'ch':
        print('Chi max: ', '{:.2e}'.format(Decimal(np.abs(CHI.min()))),)
        filename='Scalar.png'
        #Colored map of ABS(CHI)
        PLTF= np.abs(CHI)
        title = 'ABS[$\\chi$]'

    if Field == 'bs':
        print('Beta min: ', '{:.2e}'.format(Decimal(BETA.min())),)
        filename='Shift.png'
        #Colored map of BETA
        PLTF= BETA
        title = '$\\beta^\\phi$'

    if Field == 'bb':
        print('B max: ', '{:.2e}'.format(Decimal(BTOT.max()*bunit)),'[G]')
        filename='Btot.png'
        #Colored map of BTOT
        PLTF= BTOT*bunit/1.0e17
        title = '$B$' '[$10^{17}$G]'

    if Field == 'bt':
        print('Btor max: ', '{:.2e}'.format(Decimal(BTOR.max()*bunit)),'[G]')
        filename='Btor.png'
        #Colored map of BTOR
        PLTF= BTOR*bunit/1.0e17
        title = '$B_{tor}$' '[$10^{17}$G]' 

    if Field == 'bp':
        print('Bpol max: ', '{:.2e}'.format(Decimal(BPOL.max()*bunit)),'[G]')
        filename='Bpol.png'
        #Colored map of BPOL
        PLTF= BPOL*bunit/1.0e17
        title = '$B_{pol}$' '[$10^{17}$G]'
        FieldLines(NTH,NR,R,TH,Rmin,Rmax)

    if Field == 'ep':
        print('Epol/c max: ', '{:.2e}'.format(Decimal(EPOL.max()*bunit)),'[G]')
        filename='Epol.png'
        #Colored map of EPOL/c
        PLTF= EPOL*bunit/1.0e17
        title = '$E_{pol}/c$' '[$10^{17}$G]'

    if Field == 'jp':
        print('Jpol max: ', '{:.2e}'.format(Decimal(JPOL.max()*bunit)),'[G]')
        filename='Jpol.png'
        #Colored map of JPOL
        PLTF= JPOL*bunit/1.0e17/runit
        title = '$J_{pol}$' '[$10^{17}$G/km]'

    if Field == 'jt':
        print('J$_{tor}$ max: ', '{:.2e}'.format(Decimal(JTOR.max()*bunit)),'[G]')
        filename='Jtor.png'
        #Colored map of JTOR
        PLTF= JTOR*bunit/1.0e17/runit
        title = '$J_{tor}$' '[$10^{17}$G/km]'

        
    plt1=plt.pcolor(XC*runit,YC*runit,PLTF,cmap=CMAP,alpha=ALPHAC,vmin=PLTF.min(),vmax=PLTF.max())
    plt.pcolor(-XC*runit,YC*runit,PLTF,cmap=CMAP,alpha=ALPHAC,vmin=PLTF.min(),vmax=PLTF.max())
    plt.title(title,color='black',size=fs,pad=8)
    clb = plt.colorbar(plt1,fraction=fr,pad=pd,orientation="horizontal")
    clb.ax.tick_params(labelsize=lbs*0.8,colors='white')           
    clb.outline.set_edgecolor('white')
    clb.solids.set_rasterized(True)
            
    theta = np.linspace(0,np.pi,NTH)
    XCS =np.zeros(2*NTH)
    YCS =np.zeros(2*NTH)
    XCS[0:NTH] = np.sin(TH[0:NTH])*SURF[0:NTH]*runit
    YCS[0:NTH] = np.cos(TH[0:NTH])*SURF[0:NTH]*runit
    XCS[NTH:2*NTH] = -np.sin(TH[0:NTH])*SURF[0:NTH]*runit
    YCS[NTH:2*NTH] = -np.cos(TH[0:NTH])*SURF[0:NTH]*runit
    pts=(np.zeros((len(XCS),2)))
    pts[:,0]=XCS
    pts[:,1]=YCS
    tck, u = splprep(pts.T, u=None, s=0.0, per=0) 
    u_new = np.linspace(u.min(), u.max(), 1000)
    x_new, y_new = splev(u_new, tck, der=0)
    plt.plot(x_new, y_new, '-',color='white',linewidth=2)

    #Showing
    if (PLOT):
        plt.gca().set_aspect('equal', adjustable='box')
        fig.savefig(os.path.join(cwd,filename), format='png', dpi=70, bbox_inches='tight')

    plt.show()

# =============================================================================
# RUN FILE AS MAIN
# =============================================================================
if __name__ == "__main__":
    main()



