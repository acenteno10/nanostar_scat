from core_star import *

'''
Purpose:
    To model the optical properties of a Au nano star with a dielectric core.
    This code writes the shape.dat file for ddscat code. 
    


'''

R = 16 #radius of core in nm
tipR = 10 #radius of cone in nm
tipH = 20 #height of cone in nm
N=0.6666 #the dipole spacing in nm
T=4 #thickness of shell in nm
pos= nanostar_six(R,tipR,tipH,N,T)
print ('shape.dat written')
print (' Edit ddscat.par for wavelengths and dielectric material')
print ('Then you can run ddscat to obtain the extinction, absorption and scattering efficiencies for the nanostar')
print ('DDSCAT places the results in qtable')
