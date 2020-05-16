import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# 100分間の測定を想定
NSIDE = 128
nmb_pix = hp.nside2npix(NSIDE)
I_lu = np.zeros(nmb_pix)
year = 60*60*24*365
temps = 2*np.pi/year # 一日あたりの黄道座標における変化の割合
hand = 60*100# 100分間の秒数
minute = 60*15
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
    return Vect

for t in range(year):
    theta_cal = np.arccos(spin_prec(t)[1]/np.sqrt(spin_prec(t)[0]**2 + spin_prec(t)[1]**2 + spin_prec(t)[2]**2))
    #phi_cal = np.arccos(spin_prec(t)[0]/np.sqrt(spin_prec(t)[0]**2+spin_prec(t)[1]**2))
    theta.append(theta_cal)

    if spin_prec(t)[0]>=0:
        phi_cal = np.arccos(spin_prec(t)[2]/np.sqrt(spin_prec(t)[0]**2+spin_prec(t)[2]**2)) + temps*t
        if phi_cal > 2*np.pi:
          phi_cal = phi_cal - 2*np.pi

    if spin_prec(t)[0]<0:
        phi_cal = 2*np.pi - np.arccos(spin_prec(t)[2]/np.sqrt(spin_prec(t)[0]**2+spin_prec(t)[2]**2)) + temps*t
        if phi_cal > 2*np.pi:
          phi_cal = phi_cal - 2*np.pi

    phi.append(phi_cal)
    pix = hp.ang2pix(NSIDE, theta_cal, phi_cal, nest=False, lonlat=False)
    I_lu[pix] += 1

#Plot hit count maps in galactic and ecliptic coordinates
hp.mollview( I_lu,
            title="Hit count map in Ecliptic coordinates",
            unit="Hit number",
            norm="hist")
plt.show()
