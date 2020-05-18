#軌道プロット
#並列計算実装
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from joblib import Parallel, delayed
import time

start = time.time()
def create_mollweide_axes():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection="mollweide")
    return ax

def rotation(time_array):
    psi = time_array*2.*np.pi/1200
    phi = time_array*2.*np.pi/11540
    eta = time_array*2.*np.pi/60*60*24*365
    vec = np.array([np.cos(np.radians(90)),np.sin(np.radians(90)),0])
    beta = np.array([np.cos(np.radians(50)),np.sin(np.radians(50)),0])
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
    vect_1 = mat_rot1.dot(vec)

    mat_rot2 = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
    vect_2 = mat_rot2.dot(vect_1)

    mat_rot3 = np.array([[1,0,0],[0,np.cos(eta),-np.sin(eta)],[0,np.sin(eta),np.cos(eta)]])
    vect_3 = mat_rot3.dot(vect_2)
    return vect_3


NSIDE = 128
NPIX = hp.nside2npix(NSIDE)
I_lu = np.zeros(NPIX)
vec = np.array([np.cos(np.radians(95)),np.sin(np.radians(95)),0])
beta = np.array([np.cos(np.radians(45)),np.sin(np.radians(45)),0])

day = 60*60*24#1日の秒数
year = day*365
times = day
time_array = np.arange(times+1)

Vect1deg = [[],[],[]]
Vect2deg = [[], []]

psi = time_array*2.*np.pi/1200
phi = time_array*2.*np.pi/11540
eta = time_array*2.*np.pi/60*60*24*365

for i in range(times):
    Vect1deg[0].append(rotation(i)[0])
    Vect1deg[1].append(rotation(i)[1])
    Vect1deg[2].append(rotation(i)[2])

#Vect1deg = np.array(Parallel(n_jobs=-1)( [delayed(rotation)(i) for i in range(times)] )).T

#Vect1deg = np.array(Vect1deg).T

pix_d = hp.vec2pix(NSIDE,Vect1deg[0],Vect1deg[1],Vect1deg[2])
#pix_d = hp.vec2pix(NSIDE,x,y,z)
pix_d
for i in range(times):
    I_lu[pix_d[i]] += 1

hp.mollview(I_lu,title="Hit count map in Ecliptic coordinates", unit="Hit number")
hp.graticule()
elapsed_time = time.time() - start
print ("計算時間:{0}".format(elapsed_time) + "[sec]")

plt.show()
