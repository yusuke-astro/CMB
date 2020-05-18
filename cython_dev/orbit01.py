#軌道プロット
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
#import loop_func as lf

start = time.time()

def create_mollweide_axes():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection="mollweide")
    return ax


def rot_psi(time_array):
    psi = time_array*2.*np.pi/600
    phi = time_array*2.*np.pi/5760
    eta = time_array*2.*np.pi/year
    vec = np.array([np.cos(np.radians(95)),np.sin(np.radians(95)),0])
    beta = np.array([np.cos(np.radians(45)),np.sin(np.radians(45)),0])
    #vec = np.array([np.cos(np.radians(90)),np.sin(np.radians(90)),0])
    #beta = np.array([np.cos(np.radians(50)),np.sin(np.radians(50)),0])
    ux = beta[0]
    uy = beta[1]
    uz = beta[2]
    c = np.cos(psi)
    s = np.sin(psi)
    mat_rot1 = np.array([
        [ux**2*(1-c)+c , ux*uy*(1-c)-uz*s , ux*uz*(1-c)+uy*s],
        [ux*uy*(1-c)+uz*s , uy**2*(1-c)+c , uz*uy*(1-c)-ux*s],
        [ux*uz*(1-c)-uy*s , uz*uy*(1-c)+ux*s , uy**2*(1-c)+c]
        ])
    vect_1 = np.dot(mat_rot1, vec)

    mat_rot2 = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
    vect_2 = np.dot(mat_rot2, vect_1)

    mat_rot3 = np.array([[1,0,0],[0,np.cos(eta),-np.sin(eta)],[0,np.sin(eta),np.cos(eta)]])
    vect_3 = np.dot(mat_rot3,vect_2)
    return vect_3

NSIDE = 128
NPIX = hp.nside2npix(NSIDE)
I_lu = np.zeros(NPIX)
Vect1deg = [[],[],[]]
Vect2deg = [[],[]]

day = 60*60*24#1日の秒数
year = day*365
times = day

time_array = np.arange(times+1)

for i in range(times+1):
    Vect1deg[0].append(rot_psi(i)[0])
    Vect1deg[1].append(rot_psi(i)[1])
    Vect1deg[2].append(rot_psi(i)[2])

"""
for i in range(times):
    keido = np.arcsin(Vect1deg[2][i])
    ido = np.arcsin(Vect1deg[1][i]/np.cos(keido))
    Vect2deg[0].append(keido)
    Vect2deg[1].append(ido)
    pix = hp.ang2pix(NSIDE, Vect2deg[0], Vect2deg[1], nest=False, lonlat=False)
    I_lu[pix] += 1
"""

Vect1deg = np.array(Vect1deg)
pix_d = hp.vec2pix(NSIDE,Vect1deg[0],Vect1deg[1],Vect1deg[2])
pix_d
for i in range(times):
    I_lu[pix_d[i]] += 1

elapsed_time = time.time() - start
print ("計算時間:{0}".format(elapsed_time) + "[sec]")

hp.mollview(I_lu,title="Hit count map in Ecliptic coordinates", unit="Hit number")
hp.graticule()
#np.savetxt('np_savetxt04.txt', I_lu)
#a = np.loadtxt('np_savetxt04.txt')

plt.show()
