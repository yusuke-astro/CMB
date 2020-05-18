import cython
import numpy as np
import healpy as hp
cimport numpy as cnp


cdef rotation(int time_array):
  cdef:
    double psi, phi, eta
    cnp.ndarray vec, beta
    cnp.ndarray vect_1, vect_2, vect_3,
    cnp.ndarray mat_rot1, mat_rot2, mat_rot3
    double ux, uy, uz, c, s
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

def get_map(int NSIDE, int times):
  cdef:
    int i, NPIX
    list Vect1deg
    cnp.ndarray I_lu
    cnp.ndarray pix

  NPIX = hp.nside2npix(NSIDE)
  I_lu = np.zeros(NPIX)

  Vect1deg = [[],[],[]]
  for i in range(times):
    Vect1deg[0].append(rotation(i)[0])
    Vect1deg[1].append(rotation(i)[1])
    Vect1deg[2].append(rotation(i)[2])
  #orbit = np.array(Vect1deg)
  pix = hp.vec2pix(NSIDE,Vect1deg[0],Vect1deg[1],Vect1deg[2])
  for i in range(times):
    I_lu[pix[i]] += 1
  return I_lu

def orbit_vector(int times):
  cdef:
    int i
    list Vect1deg

  Vect1deg = [[],[],[]]
  for i in range(times):
    Vect1deg[0].append(rotation(i)[0])
    Vect1deg[1].append(rotation(i)[1])
    Vect1deg[2].append(rotation(i)[2])
  return Vect1deg

def spin_prec(cnp.ndarray time):
  cdef:
    double alpha, beta, omega_a, omega_b, omega_orbit
    cnp.ndarray vect_prec, vect_spin, vect_orbit, Vect, Vect_2

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
  return Vect

def get_map2(int NSIDE, cnp.ndarray time):
  cdef:
    int i, NPIX
    double alpha, beta, omega_a, omega_b, omega_orbit
    cnp.ndarray vect_prec, vect_spin, vect_orbit, Vect, Vect_2, I_lu, pix
  NPIX = hp.nside2npix(NSIDE)
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
  pix = hp.vec2pix(NSIDE,Vect[2],Vect[0],Vect[1])
  I_lu = np.zeros(NPIX)
  for i in range(len(time)):
    I_lu[pix[i]] += 1
  return I_lu


def repack(int NPIX, int times, cnp.ndarray pix):
  cdef:
    int i
    cnp.ndarray I_lu
  I_lu = np.zeros(NPIX)
  for i in range(times):
    I_lu[pix[i]] += 1
  return I_lu
