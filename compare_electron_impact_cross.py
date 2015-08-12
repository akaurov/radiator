import numpy as np
import cosmolopy
import matplotlib.pyplot as plt

from radiator.figurestyle import *
# from rstat.crossections import sigmaH
define_figure_style()
# http://arxiv.org/pdf/0807.1969.pdf


from radiator.crosssections import *

E_list = np.logspace(1, 3, 100)

data = np.loadtxt("radiator/datafiles/CCC-H-ion.dat")
plt.plot(E_list, sigmaHion(E_list/6.24e11))
plt.plot(data[:,0], data[:,1], "s")
plt.xscale('log')
plt.yscale('log')

data = np.loadtxt("radiator/datafiles/CCC-H-ex.dat")
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