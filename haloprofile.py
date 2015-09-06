__author__ = 'A'

from colossus.Cosmology import *
import numpy as np
import matplotlib.pyplot as plt

setCosmology('planck13')
cosmo = Cosmology.getCurrent()

M = 1e-12
z = 100
c = cl.HaloConcentration.concentration(M, 'vir', z, model = 'diemer15')
print np.log10(M), z, c
delta_list = np.logspace(0, 8, 33)

profile = cl.HaloDensityProfile.NFWProfile(M = 1E12, mdef = 'vir', z = z, c = c)
R = np.logspace(-5, 8, 1000000)
rho = profile.density(R)/(cl.Cosmology.AST_rho_crit_0_kpc3*cosmo.Om0*(1.+z)**3)

V = 4./3.*np.pi*R**3
V[1:] -= V[:-1]

rho2V = rho**2*V
rho2V /= rho2V.sum()

plt.plot(R, rho)
plt.xscale('log')
plt.yscale('log')

plt.hist(rho, bins=delta_list, weights=rho2V, histtype='step',)
plt.xscale('log')
plt.yscale('log')