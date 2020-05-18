#軌道プロット
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import loop_func as lf


start = time.time()

def create_mollweide_axes():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection="mollweide")
    return ax


def spin_prec(time):

    alpha = np.radians(45)
    beta = np.radians(50)
    omega_a = np.pi/30/192.348
    omega_b = 0.05*np.pi/30
    vect_prec = np.array([
        [np.cos(alpha)*np.cos(omega_a*time), -np.sin(omega_a*time), np.sin(alpha)*np.cos(omega_a*time)],
        [np.cos(alpha)*np.sin(omega_a*time), np.cos(omega_a*time), np.sin(alpha)*np.sin(omega_a*time)],
        [-np.sin(alpha), 0, np.cos(alpha)]
    ])
    vect_spin = np.array([np.sin(beta)*np.cos(omega_b*time), np.sin(beta)*np.sin(omega_b*time), np.cos(beta)])
    Vect = np.dot(vect_prec,vect_spin)

    omega_orbit = 2*np.pi/(60*60*24*365)
    vect_orbit = np.array([
        [np.cos(omega_orbit*time), 0, np.sin(omega_orbit*time)],
        [0, 1, 0],
        [-np.sin(omega_orbit*time), 0, np.cos(omega_orbit*time)]
    ])

    Vect = np.dot(vect_orbit,Vect)
    #Vect_2 = np.array(Vect[2],Vect[0],Vect[1])

    return np.array(Vect)

NSIDE = 256
NPIX = hp.nside2npix(NSIDE)
day = 60*60*24#1日の秒数
year = day*365
times = day


time_array = np.arange(times+1)


Vect1deg = spin_prec(time_array)#
#pix = hp.vec2pix(NSIDE,Vect1deg[0],Vect1deg[1],Vect1deg[2])
pix = hp.vec2pix(NSIDE,Vect1deg[2],Vect1deg[0],Vect1deg[1])
map = lf.repack(NPIX, times+1, pix)
#map = lf.get_map(NSIDE, times+1)
#map = lf.get_map2(NSIDE, time_array)

elapsed_time = time.time() - start
print ("計算時間:{0}".format(elapsed_time) + "[sec]")
hp.mollview(map, title="Hit count map in Ecliptic coordinates", unit="Hit number")
hp.graticule()

#np.savetxt('np_savetxt04.txt', I_lu)
#a = np.loadtxt('np_savetxt04.txt')

plt.show()
