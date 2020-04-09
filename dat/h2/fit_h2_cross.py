"""script that was used to come up with the approximation for the H2 ionization cross section at < 17 nm"""

import numpy as np
import matplotlib.pyplot as plt

PLANCK = 6.62606957e-27
LIGHT = 2.99792458e10
NM = 1.0e-7
EV = 1.60218e-12


def h_cs(w):
    """ionization cross section of H, ported from Ionization.cpp. 

    w: wavelength in nanometer

    """
    E0 = 4.298e-1
    sigma0 = 5.475e4
    ya = 3.288e1
    P = 2.963
    frequency = LIGHT / (w * NM)
    x = PLANCK * frequency / (E0 * EV)
    y = x

    Fy = (
        (x - 1)
        * (x - 1)
        * np.power(y, 0.5 * P - 5.5)
        * np.power(1 + np.sqrt(y / ya), -P)
    )

    return sigma0 * Fy * 1e-18


wav, ion_cs = np.loadtxt("crossections_leiden.txt", usecols=(0, 3), unpack=True)

fig, axs = plt.subplots(2, 1, sharex=True)

h2 = ion_cs
h = h_cs(wav)
ratio = h2 / h

# now fit something using data < 70 nm
lambda_max = 28
fit_where = np.where(wav < lambda_max)

# power law?
# fit_x = np.log(wav[fit_where])
# fit_y = np.log(h2[fit_where])

# simplest fit: average of the ratio
avg_ratio = np.average(ratio[fit_where])
print(
    "average ratio h2 / h cross section = {} below {} nm".format(avg_ratio, lambda_max)
)

axs[0].loglog(wav, ion_cs, label="h2")
axs[0].loglog(wav, h_cs(wav), label="h")
axs[0].loglog(wav, avg_ratio * h_cs(wav), label="h times ratio")
axs[1].semilogx(wav, ratio)
axs[1].axhline(avg_ratio)
plt.show()
