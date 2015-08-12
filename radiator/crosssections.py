__author__ = 'kaurov'
import numpy as np


sigmaT = 6.6524e-25 #cm^2
me = 9.109e-28 # g
c = 2.9979e10 # cm/s
#hbar = 4.135668e-15/2.0/np.pi # eV*s
hbar = 1.0545e-27 # erg*s
kB = 1.380648e-16 # erg/K

#
# Photoionization cross sections
#

# http://adsabs.harvard.edu/abs/1996ApJ...465..487V

lookup = np.zeros([10, 10, 20])
lookup[1, 1, :9] = [1.360e1,    5.000E+4,   4.298E-1, 5.475e+4, 3.288e1,    2.963e0,    0.000e0,    0,          0]
lookup[2, 2, :9] = [24.59,      5e4,        13.61,    949.2,    1.469,      3.188,      2.039,      0.4434,     2.136]
lookup[2, 1, :9] = [54.42,      5e4,        1.72,     13690,    32.88,      2.963,      0,          0,          0]
lookup[8, 8, :9] = [1.362e1,    5.380E+2,   1.240e0,  1.745E+3, 3.784e0,    1.764e1,    7.589E-2,   8.698e0,    1.271E-1]

def sigmaX(E, Z, N):
    '''
    Photoionization cross section
    :param E: energy of a photon
    :param Z: Atomic number
    :param N: number of electrons
    :return: cross section in cm^-2
    '''
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


#
# Electron impact cross sections
#

def sigmaHex(E):
    '''
    Excitation cross section
    :param E: Energy of incident electron in ergs
    :return: cm^-2
    '''
    return (E*6.24e11 > 10.2) * 3.75e-15 * np.log(E/2.*6.24e11) / (E*6.24e11)

# Shull 1985
def sigmaHion(E):
    '''
    Ionization cross section
    :param E: Energy of incident electron in ergs
    :return: cm^-2
    '''
    return 2.75e-15 * np.log(E/13.6*6.24e11) / (E*6.24e11)

# CCC fits
def sigmaHe(E, type = 'ion', ion_state='I'):
    '''
    Ionization and excitation cross sections for HeI and HeII
    :param E: Energy of incident electron in ergs
    :return: cm^-2
    '''
    if (type == 'ion') and (ion_state == 'I'):
        return 2.75e-15 * np.log(E/24.6*6.24e11) / (E*6.24e11)
    elif (type == 'ion') and (ion_state == 'II'):
        return 0.75e-15 * np.log(E/54.4*6.24e11) / (E*6.24e11)
    elif (type == 'ex') and (ion_state == 'I'):
        return (E*6.24e11 > 20) * 3.75e-15 * np.log(E/3.*6.24e11) / (E*6.24e11)
    elif (type == 'ex') and (ion_state == 'II'):
        return (E*6.24e11 > 41) * 0.025e-15 * np.log(E/4.6*6.24e11) / (E*6.24e11)
    else:
        return "error"

# Secondary electron energy distribution
def rhoE(E, ei):
    '''
    Energy of secondary electron
    :param E: Energy of primary electron in ergs
    :param ei: ei in ergs
    :return: Energy of secondary electron in ergs
    '''
    r = np.random.rand(len(E))
    ei_mod = ei+np.log10(E*6.24e11/100)*2/6.24e11
    temp = np.tan(r*np.pi/2.0)*(ei*(1+np.log10(E*6.24e11/ei_mod)))
    temp[temp > E/2.0] = E/2.0
    return temp

# Electron cooling due to interaction with other electrons
# TODO Add a better Coulomb logarithm fitting function
def eedEdt(E, ne, T):
    '''
    Rate of energy loss of an electron with energy E in electron plasma with density ne and temperature T
    :param E: in ergs
    :param ne: in cm^-3
    :param T: in K
    :return: in ergs/s
    '''
    omega = c * np.sqrt(1.0 - 511e3**2 / (511e3 + E*6.24e11)**2)
    lnL = 16.3 - np.log10(ne)*1.15 + 3.45*np.log10(T/100.0) # my fit to Spitzer 1965
    return 4. * np.pi * ne * 4.8e-10**4 / 9.1e-28 / omega * lnL

def sigmakn(Eg, e, gamma):
    Gamma = 4*e*gamma/me/c**2
    eta = e*Eg/(me*c**2)**2
    q = Eg/Gamma/(gamma*me*c**2-Eg)
    G = 2.0*q*np.log(q)+(1.0+2.0*q)*(1.0-q)+2.0*eta*q*(1.0-q)
    G[G < 0] = 0
    G[(4*gamma**2)**-1 > q] = 0
    G[q > 1] = 0
    return 3.0*sigmaT/4.0/e/gamma**2*G