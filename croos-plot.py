__author__ = 'kaurov'
import numpy as np
import matplotlib.pyplot as plt
from radiator.figurestyle import *
# from rstat.crossections import sigmaH
define_figure_style()
# http://arxiv.org/pdf/0807.1969.pdf



from radiator.crosssections import *

sigmaT = 6.6524e-25 #cm^2
me = 9.109e-28 # g
c = 2.9979e10 # cm/s
#hbar = 4.135668e-15/2.0/np.pi # eV*s
hbar = 1.0545e-27 # erg*s
kB = 1.380648e-16 # erg/K


Eg_list = np.logspace(-2, 12, 1000)/6.24e11
temp = np.zeros([len(Eg_list), 4])
EgeV_list = Eg_list*6.24e11

from scipy.interpolate import InterpolatedUnivariateSpline

crossections = np.loadtxt('radiator/datafiles/crossdata.dat')
# crossectionsHe = np.loadtxt('radiator/datafiles/crossdataHe.dat')
# crossectionsO = np.loadtxt('radiator/datafiles/crossdataO.dat')
# crossectionsHe2 = np.loadtxt('radiator/datafiles/photoionHe.dat')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-3]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, ':', lw=2, label='Electron pr. prod.')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-4]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, ':', lw=2, label='Nucl. pr. prod.')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-7]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '-.', lw=2, label='scattering coher.')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-6]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '--', lw=2, label='scattering incoher.')
# s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
# factor = 10**s(np.log10(EgeV_list/1e6))
# plt.plot(EgeV_list, factor*1e-24, '-', lw=1, label='O cross total')
plt.plot(EgeV_list, sigmaX(EgeV_list, 1, 1)*1e-18, lw=2, label='Photoel.')
# plt.plot(EgeV_list, sigmaX(EgeV_list, 8, 8)*1e-18, label='O phot 1996')
plt.xscale('log')
plt.xlabel(r'$E\;\mathrm{[eV]}$')
plt.ylabel(r'$\sigma\;\mathrm{[cm^{2}]}$')
plt.yscale('log')
plt.legend()