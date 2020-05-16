import matplotlib.pyplot as plt
import numpy as np
import math
import healpy as hp
from scipy.optimize import curve_fit
from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D

def create_mollweide_axes():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection="mollweide")
    return ax

def rot_x(Vect,phi):
    mat_rot1 = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
    vect_2 = mat_rot1.dot(Vect)
    return vect_2

def rot_psi(Vect,U,psi):
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

tmps = 60*60*24#1日の秒数
year_sec = tmps*365
degree = 360



Vect = [np.cos(np.radians(95)),np.sin(np.radians(95)),0]
U = [np.cos(np.radians(45)),np.sin(np.radians(45)),0]
rot_psi(Vect,U,np.pi)


psi = np.linspace(0, 2*np.pi, degree)
psi

phi = []
psi = []
creve = []


for t in range(tmps):
    psi.append(2*np.pi*t/600)
    phi.append(2*np.pi*t/5760)
    creve.append(2*np.pi*t/year_sec)
Vect1deg=[[],[],[]]

for i in range(tmps):
    rott1 = rot_psi(Vect,U,psi[i])
    rott2 = rot_x(rott1,phi[i])
    Vect1deg[0].append(rott2[0])
    Vect1deg[1].append(rott2[1])
    Vect1deg[2].append(rott2[2])
Vect2deg = [[], []]

r = 1.5e6
for i in range(tmps):
    beta = np.arcsin(Vect1deg[2][i])
    lambd = np.arcsin(Vect1deg[1][i]/np.cos(beta))
    Vect2deg[0].append(beta)
    Vect2deg[1].append(lambd)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Vect1deg[0],Vect1deg[1],Vect1deg[2])

#plt.figure()
create_mollweide_axes()
plt.plot(Vect2deg[0], Vect2deg[1], '-', color="red", alpha=0.5)
plt.title("Trajectory of the Planck spin axis for 1 year in Galactic Coordinates");

plt.show()
