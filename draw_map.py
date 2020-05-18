import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

print("Input file path...")
print(">>>")
file_path = input()
I_lu = np.loadtxt(file_path)
I_lu
hp.mollview(I_lu, title="Hit count map in Ecliptic coordinates", unit="Hit number", norm="hist",cmap="jet")
plt.show()
