__author__ = 'A'

from colossus.Cosmology import *
from colossus.HaloConcentration import *
from colossus.HaloDensityProfile import *
import numpy as np
import matplotlib.pyplot as plt

setCosmology('planck13')
cosmo = getCurrent()

M = 1e-6
z = 100
c = concentration(M, 'vir', z, model='diemer15')
print np.log10(M), z, c
delta_list = np.logspace(0, 8, 33)

profile = NFWProfile(M=M, mdef='vir', z=z, c=c)
R = np.logspace(-12, 8, 1000000)
rho = profile.density(R)/(AST_rho_crit_0_kpc3*cosmo.Om0*(1.+z)**3)

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
plt.ylim([1e-5,1])