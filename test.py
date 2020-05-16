#軌道プロット

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def rot_x(Vect,phi):
    mat_rot1 = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
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

vec = np.array([np.cos(np.radians(95)),np.sin(np.radians(95)),0])
beta = np.array([np.cos(np.radians(45)),np.sin(np.radians(45)),0])

"""
N = 50
psi = np.linspace(0, 2*np.pi,N)
Vect1deg = [[],[],[]]

for i in range(0,N):
    roted = rot_psi(vec,beta,psi[i])

    Vect1deg[0].append(roted[0])
    Vect1deg[1].append(roted[1])
    Vect1deg[2].append(roted[2])
"""
day_sec = 60*60*24#1日の秒数
year_sec = day_sec*365
phi = []
psi = []
creve = []
Vect1deg=[[],[],[]]

for t in range(day_sec):
    psi.append(2*np.pi*t/600)#600=10分で１週
    phi.append(2*np.pi*t/5760)#5760=96分でちk
    #creve.append(2*np.pi*t/year_sec)

for i in range(day_sec):
    rott1 = rot_psi(vec,beta,psi[i])
    rott2 = rot_x(rott1,phi[i])
    Vect1deg[0].append(rott2[0])
    Vect1deg[1].append(rott2[1])
    Vect1deg[2].append(rott2[2])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Vect1deg[0],Vect1deg[1],Vect1deg[2],"-")



"""
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(0, 0, 0, 1,0,0, length=0.1, normalize=True,color="red")
ax.text(0.5, 0, 0, "x", color = "red", size = 15)

ax.quiver(0, 0, 0, 0,1,0,length=0.1, normalize=True,color="red")
ax.text(0, 0.5, 0, "y", color = "red", size = 15)

ax.quiver(0, 0, 0, 0,0,1, length=0.1, normalize=True,color="red")
ax.text(0, 0, 0.5, "z", color = "red", size = 15)

ax.quiver(0, 0, 0, vec[0],vec[1],vec[2], length=0.1, normalize=True)
ax.quiver(0, 0, 0, vec1[0],vec1[1],vec1[2], length=0.1, normalize=True, color="green")
ax.quiver(0, 0, 0, beta[0],beta[1],beta[2], length=0.1, normalize=True, color="orange")

ax.set_xlabel("x", fontsize = 16)
ax.set_ylabel("y", fontsize = 16)
ax.set_zlabel("z", fontsize = 16)
"""
#plt.show()
