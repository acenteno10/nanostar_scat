# This python script creates a shape.dat file for a 6 tip nanostar #
import numpy as np
import math
from numpy import random
from scipy.spatial.transform import Rotation



#This generates a 6 spike nanostar with a core.
def nanostar_six(R,coneR,coneH,N,T):
    pos = core(R,T,N)
    pos = unique_rows(pos)
# first cone
    con=[]    
    con = cone(coneR,coneH,N,(R+T))
#other cones
    con1=[]
    con2=[]
    con3=[]
    con4=[]
    con5=[]
# Euler rotation to generate 6 cones.
    r = Rotation.from_euler('y', 90, degrees=True)
    con1 = r.apply(con)
    con2=r.apply(con1)
    con3=r.apply(con2)
    r = Rotation.from_euler('x', 90, degrees=True)
    con4=r.apply(con)
    r = Rotation.from_euler('x', 270, degrees=True)
    con5=r.apply(con)   
    pos = np.concatenate((pos,con,con1,con2,con3,con4,con5))
    pos = np.array(pos).astype(int)
    pos = list(pos)
    pos = unique_rows(pos)
    writeDDSCAT(pos,R,T,N,coneH)
    return(pos)

def cone(coneR,coneH,N,R=0,T=0):
    coneT = np.arctan(coneR/coneH)
    hyp = np.sqrt(coneR**2+coneH**2)
    coneR = np.sin(coneT)*hyp
    T1 = 0 #T1 is the depth of cone into the sphere so the end points sit on sphere surface
    R=R+T #Make sure that r is radius of sphere including shell and core
    if R != 0:
        T1 = -np.sqrt(R**2-coneR**2)+R
    H = coneH+T1
    cone = []
    for xi in np.linspace(int(-coneR/N), int(coneR/N),int(coneR*2/N)+1):
        for yi in np.linspace(int(-coneR/N), int(coneR/N),int(coneR*2/N)+1):
#          print (H,N)
          for zi in np.linspace(int(0),int(H/N),int(H/N)+1):
                if ((abs(zi) < H/N) and (np.sqrt(xi**2+yi**2) < (H/N-zi)/(float(H/N))*(coneR/N))):
                    cone.append([xi,yi,zi+int(R/N-(T1/N))])
    return np.array(cone)

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

#This bit generates the shape.dat file for running DDSCAT and writes aeff to ddscat.par
def writeDDSCAT(pos,R,T,N,coneH=0):
    Ls=(R+T+coneH)
    Li=0
    if T>0: Li=(R)    
    tot = pos
    nParticles = len(tot)
    x = int(np.max(tot[0,:]))- int(np.min(tot[0,:]))
    y = int(np.max(tot[1,:]))- int(np.min(tot[1,:]))
    z = int(np.max(tot[2,:]))- int(np.min(tot[2,:]))
    with open( 'shape.dat', 'w' ) as g:
        g.write(' >TARREC   rectangular prism; AX,AY,AZ= ' +repr(x)+' '+repr(y)+' '+repr(z)+'\n     '+repr(nParticles)+' = NAT \n')
        g.write('  1.000000  0.000000  0.000000 = A_1 vector\n')
        g.write('  0.000000  1.000000  0.000000 = A_2 vector\n')
        g.write('  1.000000  1.000000  1.000000 = lattice spacings (d_x,d_y,d_z)/d\n')
        g.write('  0.000000  0.000000  0.000000 = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d for dipole 0 0 0\n')
        g.write('       JA  IX  IY  IZ ICOMP(x,y,z)\n')
        for t in range(len(pos)):
            xi = pos[t,0]
            yi = pos[t,1]
            zi = pos[t,2]
            rad=(np.sqrt(xi**2+yi**2+zi**2))
            if rad <= Ls/N and rad > Li/N: g.write( "      {0:.0f}   {1:.0f}   {2:.0f}   {3:.0f} {4:.0f} {5:.0f} {6:.0f}\n".format( t+1, pos[ t, 0 ], pos[ t, 1 ], pos[ t, 2 ], 2, 2, 2 ) )
            elif (rad <= Li/N): g.write( "      {0:.0f}   {1:.0f}   {2:.0f}   {3:.0f} {4:.0f} {5:.0f} {6:.0f}\n".format( t+1, pos[ t, 0 ], pos[ t, 1 ], pos[ t, 2 ], 1, 1, 1 ) )
    v = len(pos)
    aeff = (v*(3/4)/np.pi)**(1/3)*N
    with open('_ddscat.par', 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('REPRAD',"{:.4f}".format(aeff/1000))
    with open('ddscat.par', 'w') as file:
        file.write(filedata)
    return


def core(R,T,N):
#	This code defines the sphere with a dielectric core (refractive index of core >1). R is radius of core and T is 
#       thickness of the shell
	ref=[]
	Ls=(R+T)
	Li=(R)
#	print (2*Ls/N)
	for xi in np.linspace(int (-Ls/N), int (Ls/N), int ((2*Ls/N)+1)):
#		print (xi)
		for yi in np.linspace(int(-Ls/N), int (Ls/N),int(2*Ls/N)+1):
			for zi in np.linspace(int (-Ls/N), int (Ls/N),int(2*Ls/N)+1):
				rad=(np.sqrt(xi**2+yi**2+zi**2))
				if rad <= Ls/N: ref.append([xi,yi,zi])
	
	ref = np.array(ref).astype(int)
	pos = unique_rows(ref)
	return(pos)
