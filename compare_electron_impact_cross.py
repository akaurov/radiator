import numpy as np
import cosmolopy
import matplotlib.pyplot as plt

from radiator.figurestyle import *
# from rstat.crossections import sigmaH
define_figure_style()
# http://arxiv.org/pdf/0807.1969.pdf


from radiator.crosssections import *

E_list = np.logspace(1, 9, 1000)

data = np.loadtxt("radiator/datafiles/CCC-H-ion.dat")
plt.plot(E_list, sigmaHion(E_list/6.24e11),':', label="Schull '85")
plt.plot(E_list, sigmaBEQ(E_list, 13.6057, 13.6057, 1), '--', label="$\sigma_{\mathrm{BEQ}}$")
plt.plot(E_list, sigmaRBEQ(E_list, 13.6057, 13.6057, 1), '-', lw=1, label="$\sigma_{\mathrm{RBEQ}}$")
plt.plot(E_list, sigma_AR(E_list, mode='HI'), '-', lw=1, label="$AR$")
plt.plot(data[:,0], data[:,1], "s")
plt.legend(loc=3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E\;\mathrm{[eV]}$')
plt.ylabel('$\sigma\;\mathrm{[cm^{-2}]}$')

data = np.loadtxt("radiator/datafiles/CCC-H-ex.dat")
plt.plot(E_list, sigmabe(E_list, 10.2, 13.6, 1), '--', label="$\sigma_{\mathrm{BEQ}}$")
plt.plot(E_list, sigmaHex(E_list/6.24e11))
plt.plot(data[:,0], data[:,1], "s")
plt.xscale('log')
plt.yscale('log')

data = np.loadtxt("radiator/datafiles/CCC-HeI-ion.dat")
plt.plot(E_list, sigmaHe(E_list/6.24e11, type='ion', ion_state='I'))
plt.plot(data[:, 0], data[:, 1], "s")
plt.xscale('log')
plt.yscale('log')

data = np.loadtxt("radiator/datafiles/CCC-HeI-ex.dat")
plt.plot(E_list, sigmaHe(E_list/6.24e11, type='ex', ion_state='I'))
plt.plot(data[:, 0], data[:, 1], "s")
plt.xscale('log')
plt.yscale('log')

data = np.loadtxt("radiator/datafiles/CCC-HeII-ion.dat")
plt.plot(E_list, sigmaHe(E_list/6.24e11, type='ion', ion_state='II'))
plt.plot(data[:, 0], data[:, 1], "s")
plt.xscale('log')
plt.yscale('log')

data = np.loadtxt("radiator/datafiles/CCC-HeII-ex.dat")
plt.plot(E_list, sigmaHe(E_list/6.24e11, type='ex', ion_state='II'))
plt.plot(data[:, 0], data[:, 1], "s")
plt.xscale('log')
plt.yscale('log')