#軌道プロット

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def create_mollweide_axes():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection="mollweide")
    return ax

def rot_x(Vect,phi):
    mat_rot1 = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
    vect_2 = mat_rot1.dot(Vect)
    return vect_2

def rot_z(Vect,phi):
    mat_rot1 = np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])
    vect_2 = mat_rot1.dot(Vect)
    return vect_2

def rot_psi(Vect,U,psi):#U=beta
    ux = U[0]
    uy = U[1]
    uz = U[2]
    c = np.cos(psi)
    s = np.sin(psi)
    mat_rot2 = np.array([
        [ux**2*(1-c)+c , ux*uy*(1-c)-uz*s , ux*uz*(1-c)+uy*s],
        [ux*uy*(1-c)+uz*s , uy**2*(1-c)+c , uz*uy*(1-c)-ux*s],
        [ux*uz*(1-c)-uy*s , uz*uy*(1-c)+ux*s , uy**2*(1-c)+c]
        ])
    vect_2 = mat_rot2.dot(Vect)
    return vect_2

NSIDE = 128
NPIX = hp.nside2npix(NSIDE)
I_lu = np.zeros(NPIX)
vec = np.array([np.cos(np.radians(95)),np.sin(np.radians(95)),0])
beta = np.array([np.cos(np.radians(45)),np.sin(np.radians(45)),0])

day = 60*60*24#1日の秒数
year = day*365
times = day
phi = []
psi = []
eta = []

Vect1deg = [[],[],[]]
Vect2deg = [[], []]
map_2D = [[],[]]
for t in range(times):
    psi.append(2*np.pi*t/600)#600=10分で１週
    phi.append(2*np.pi*t/5760)#5760=96分でちk
    eta.append(2*np.pi*t/year)
    #creve.append(2*np.pi*t/times)

for i in range(times):
    rott1 = rot_psi(vec,beta,psi[i])
    rott2 = rot_x(rott1,phi[i])
    rott3 = rot_z(rott2,eta[i])
    Vect1deg[0].append(rott3[0])
    Vect1deg[1].append(rott3[1])
    Vect1deg[2].append(rott3[2])

"""
for i in range(times):
    keido = np.arcsin(Vect1deg[2][i])
    ido = np.arcsin(Vect1deg[1][i]/np.cos(keido))
    Vect2deg[0].append(keido)
    Vect2deg[1].append(ido)
    #pix = hp.ang2pix(NSIDE, Vect2deg[0], Vect2deg[1], nest=False, lonlat=False)
    #I_lu[pix] += 1
"""
#max(I_lu)
#Vect1deg[2]
Vect1deg = np.array(Vect1deg)
#ANG = hp.vec2ang(Vect1deg)
#pix = hp.ang2pix(NSIDE,ANG[0],ANG[1])
pix_d = hp.vec2pix(NSIDE,Vect1deg[0],Vect1deg[1],Vect1deg[2])
for i in range(times):
    I_lu[pix_d[i]] += 1

plt.figure()
hp.mollview(I_lu,title="Hit count map in Ecliptic coordinates", unit="Hit number")
hp.graticule()
np.savetxt('np_savetxt04.txt', I_lu)
a = np.loadtxt('np_savetxt04.txt')

"""
Theta = np.array(Vect2deg[0])
Phi = np.array(Vect2deg[1])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Vect1deg[0],Vect1deg[1],Vect1deg[2],"-")

create_mollweide_axes()
plt.plot(ANG[0],ANG[1], '-', color="red", alpha=0.5)
plt.title("Trajectory of the Planck spin axis for 1 year in Galactic Coordinates")
"""
plt.show()
