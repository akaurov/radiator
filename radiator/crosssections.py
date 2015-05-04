__author__ = 'kaurov'
import numpy as np

# http://adsabs.harvard.edu/abs/1996ApJ...465..487V

lookup = np.zeros([10, 10, 20])
lookup[1, 1, :9] = [1.360e1,    5.000E+4,   4.298E-1, 5.475e+4, 3.288e1,    2.963e0,    0.000e0,    0,          0]
lookup[2, 2, :9] = [24.59,      5e4,        13.61,    949.2,    1.469,      3.188,      2.039,      0.4434,     2.136]
lookup[2, 1, :9] = [54.42,      5e4,        1.72,     13690,    32.88,      2.963,      0,          0,          0]
lookup[8, 8, :9] = [1.362e1,    5.380E+2,   1.240e0,  1.745E+3, 3.784e0,    1.764e1,    7.589E-2,   8.698e0,    1.271E-1]

def sigmaX(E, Z, N):
    Eth, Emax, E0, sigma0, ya, P, yw, y0, y1 = lookup[Z, N, :9]
    x = E/E0-y0
    y = np.sqrt(x**2+y1**2)
    F = ((x-1.0)**2+yw**2)*y**(0.5*P-5.5)*(1.0+np.sqrt(y/ya))**(-P)
    F[E<Eth] = 0.0
    return sigma0*F

# plt.plot(EgeV_list, sigmaX(EgeV_list, 8, 8)/sigmaX(EgeV_list, 1, 1), '-', label=r'$\mathrm{OI}$'); plt.yscale('log'); plt.xscale('log')
# plt.plot(EgeV_list, sigmaX(EgeV_list, 2, 2)/sigmaX(EgeV_list, 1, 1), '--', label=r'$\mathrm{HeI}$'); plt.yscale('log'); plt.xscale('log')
# plt.plot(EgeV_list, sigmaX(EgeV_list, 2, 1)/sigmaX(EgeV_list, 1, 1), '-.', label=r'$\mathrm{HeII}$'); plt.yscale('log'); plt.xscale('log')
# plt.legend(loc=4)
# plt.xlabel(r'$\mathrm{eV}$')
# plt.ylabel(r'$\sigma_X/\sigma_H$')
# fixlogax(plt.gca())
# fixlogax(plt.gca(), 'y')
# plt.xlim([10,1e5])
# plt.ylim([0.1,1e3])