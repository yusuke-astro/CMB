import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import time as TIME

start = TIME.time()
# 3年間の測定を想定
NSIDE = 256
nmb_pix = hp.nside2npix(NSIDE)
I_lu = np.zeros(nmb_pix)
hand = 60*100# 100分間の秒数
minute = 60
year = 60*60*24*365
day = 60*60*24
month = 60*60*24*30
theta = []
phi = []

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

    return Vect

time = year*3 #この変数に計算したい測定期間をいれる
time_2 = np.arange(0,time+1,1)

vec = spin_prec(time_2)
vec
pix = hp.vec2pix(NSIDE, vec[2], vec[0], vec[1])


for i in range(time):
    I_lu[pix[i]] += 1

elapsed_time = TIME.time() - start
print ("計算時間:{0}".format(elapsed_time) + "[sec]")
hp.mollview( I_lu,
            title="Hit count map in Ecliptic coordinates",
            unit="Hit number",
            norm="hist")


plt.show()
