
def CLfit(z):
    return 1.0+8e5*z**-1.4*np.exp(-(z/40)**2)

import numpy as np
import matplotlib.pyplot as plt
from radiator.figurestyle import *
from radiator.crosssections import *
from colossus import Cosmology

sigmaT = 6.6524e-25 #cm^2
me = 9.109e-28 # g
c = 2.9979e10 # cm/s
#hbar = 4.135668e-15/2.0/np.pi # eV*s
hbar = 1.0545e-27 # erg*s
kB = 1.380648e-16 # erg/K

cs = cross_sections(cs={'photion': 'VFKY1996', 'collion': 'AR', 'collex': 'SKD'})

from scipy.interpolate import InterpolatedUnivariateSpline

crossections = np.loadtxt('radiator/datafiles/crossdata.dat')
crossectionsHe = np.loadtxt('radiator/datafiles/crossdataHe.dat')
crossectionsO = np.loadtxt('radiator/datafiles/crossdataO.dat')
crossectionsHe2 = np.loadtxt('radiator/datafiles/photoionHe.dat')


Eg_list = np.logspace(-2, 12, 1000)/6.24e11
temp = np.zeros([len(Eg_list), 4])
EgeV_list = Eg_list*6.24e11

dm_data = np.genfromtxt('PPPC4DM/AtProductionNoEW_positrons.dat', names=True)
def get_electron_spec(e=5):
    '''
    Returns the electron energy specrum
    '''
    return 10**dm_data['Log10x'][dm_data['mDM']==e], dm_data['q'][dm_data['mDM']==e]

def ICon30CMB(EgeV_list, E0, z):
    return (EgeV_list/1e9)**0.5*5e-2 *(np.exp(-(EgeV_list/1e9*3.4e3*(1e9/E0/((1.+z)/61)**0.5)**2)))


cosmo = Cosmology.setCosmology('planck15')

# z_list_ann = z_list.copy()
# boost_ann = 1+rho_dm_2[:]
# boost_ann = CLfit(z_list)

# Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11
z_start = 10000.0
z_step = 0.01
z_list = np.array([0])
f = (Eg_list[1]/Eg_list[0])**1
double_shift = int(np.ceil(np.log(2.0)/np.log(f)))
for i in range(1000):
    z_list = np.append(z_list, f*(1.+z_list[-1])-1.0)
z_list = z_list[z_list <= z_start]
z_list = z_list[::-1]
# boost = np.interp(z_list, z_list_ann, boost_ann)
boost = CLfit(z_list)

import cosmolopy as cp
age_list = cp.distance.age(z_list,**cp.fidcosmo)

electrons = np.zeros(len(EgeV_list))
electrons = EgeV_list / 6.24e11

DM_e=5
x, y = get_electron_spec(e=DM_e)
p2 = np.interp(EgeV_list, x*DM_e*1e9, y)
phot2 = 0#np.interp(EgeV_list, phot[:,0]*1e9, phot[:,1], left=0, right=0)

z = z_list[0]
T_CMB = 2.7*(1+z)
CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

tau = 3.14e7*1e2
temp = 0

Sigma_pair_prod = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,4]+crossections[:,5]))*1e-24
Sigma_pair_prod[np.isnan(Sigma_pair_prod)] = 0.0

s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-1]), k=1)
Sigma_photoion = 10**s(np.log10(EgeV_list/1e6))*1e-24
Sigma_photoion[np.isnan(Sigma_photoion)] = 0.0

Sigma_collisional_ion = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,2]))*1e-24
Sigma_collisional_ion[np.isnan(Sigma_collisional_ion)] = 0.0

s = InterpolatedUnivariateSpline(np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,-1]), k=1)
Sigma_photoion_He = 10**s(np.log10(EgeV_list/1e6))*1e-24
Sigma_photoion_He[np.isnan(Sigma_photoion_He)] = 0.0

# Sigma_photoion_O = sigmaX(EgeV_list, 8, 8)*1e-18

Sigma_collisional_ion_He = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,2]))*1e-24
Sigma_collisional_ion_He[np.isnan(Sigma_collisional_ion_He)] = 0.0

Sigma_collisional_scatter = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,1]))*1e-24
Sigma_collisional_scatter[np.isnan(Sigma_collisional_scatter)] = 0.0

# s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
# Sigma_photoion_O = 10**s(np.log10(EgeV_list/1e6))*1e-24
# Sigma_photoion_O[np.isnan(Sigma_photoion_He)] = 0.0

n_ion = np.zeros([len(z_list), 10])
electronsN = np.zeros(len(EgeV_list))
photons_particles = np.zeros(len(Eg_list))

for ii in range(len(z_list)-1):
    # plt.clf()
    # # plt.subplot(211)
    # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2,'-')
    # plt.plot(EgeV_list, electronsN/np.gradient(Eg_list)*Eg_list**2,'--')
    # plt.ylim([1e15,1e22])
    # plt.subplots_adjust(top=0.9)
    # plt.xscale('log')
    # plt.xlabel('eV')
    # plt.yscale('log')
    print ii
    dt = age_list[ii+1]-age_list[ii]
    tau = dt
    print tau
    z = 0.5 * (z_list[ii+1] + z_list[ii])
    plt.title('z=%5.1f'%z)
    T_CMB = 2.7*(1+z)
    # redshift the system
    photons_particles[:-1] = photons_particles[1:]
    photons_particles[-1] = 0.
    # update CMB
    CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
    # add sources
    photons_particles += (phot2*np.gradient(np.log10(EgeV_list)))*boost[ii]*(1.0+z)**3*tau
    electronsN = 2*(p2*np.gradient(np.log10(EgeV_list)))*boost[ii]*(1.0+z)**3*tau
    n_ion[ii,9] = (1.0+z)**3*tau*boost[ii]
    # add IC photons
    electronsN[EgeV_list<13.6] = 0
    n_ion[ii, 0] += np.sum(electronsN[EgeV_list<1e6]*EgeV_list[EgeV_list<1e6])/13.6*0.3
    electronsN[EgeV_list<1e6] = 0
    for i in range(len(electrons)):
        photons_particles += electronsN[i] * ICon30CMB(EgeV_list, electrons[i]*6.24e11, z) * np.gradient(Eg_list) / Eg_list**2
    electronsN *= 0
    # column density
    N_naked = 2.2e-7*c*tau*(1+z)**3
    N_H = N_naked*0.76
    N_He = N_naked*0.23/4
    N_O = N_naked*0.01/8
    # pair production
    temp_factor = (Sigma_pair_prod*N_naked)
    temp_factor[temp_factor > 1.0] = 1.0
    temp = photons_particles * temp_factor
    photons_particles -= temp
    electronsN[:-double_shift] += temp[double_shift:]*2
        # photoion
    temp_factor_H = (Sigma_photoion*N_H)
    temp_factor_He = (Sigma_photoion_He*N_He)
    # temp_factor_O = (Sigma_photoion_O*N_O)
    temp_factor_H[EgeV_list < 13.6] = 0.0
    temp_factor_He[EgeV_list < 24.6] = 0.0
    # temp_factor_O[EgeV_list < 13.6] = 0.0
    temp_factor = temp_factor_H + temp_factor_He# + temp_factor_O
    temp_factor_H[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor_He[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    # temp_factor_O[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor[temp_factor > 1.0] = 1.0
    temp = photons_particles * temp_factor
    n_ion[ii, 0] += np.sum(photons_particles * temp_factor_H)
    n_ion[ii, 4] += np.sum(photons_particles * temp_factor_He)
    # n_ion[ii, 8] += np.sum(photons_particles * temp_factor_O)
    photons_particles -= temp
    # plt.subplot(212)
    # plt.plot(EgeV_list, temp_factor_H, '-', lw=2)
    # plt.plot(EgeV_list, temp_factor_He, '-', lw=1)
    # plt.xscale('log')
    # plt.yscale('log')
    electronsN += temp
    # collisional ion
    temp_factor_H = (Sigma_collisional_ion*N_H)
    temp_factor_He = (Sigma_collisional_ion*N_He)
    temp_factor_H[EgeV_list < 13.6] = 0.0
    temp_factor_He[EgeV_list < 24.6] = 0.0
    temp_factor = temp_factor_H + temp_factor_He
    temp_factor_H[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor_He[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor[temp_factor > 1.0] = 1.0
    temp = photons_particles * temp_factor
    photons_particles -= temp
    n_ion[ii, 1] += np.sum(photons_particles * temp_factor_H)
    n_ion[ii, 5] += np.sum(photons_particles * temp_factor_He)
    # bad approx!!!!
    electronsN += temp
    # plt.savefig('figures/DMspectraani/test%04d.png'%ii, dpi=250)

# plt.plot(z_list, n_ion[:, 0])
# plt.plot(z_list, n_ion[:, 4])
# plt.plot(z_list, n_ion[:, 8])
# plt.xscale('log')
# plt.yscale('log')

sigmav = 2.16e-26
mx = DM_e*1e9
nb = 2.2e-7

Factor = ((Cosmology.AST_rho_crit_0_kpc3*cosmo.Om0 / Cosmology.AST_kpc**3 * Cosmology.AST_Msun) / (mx * 1.783e-33)) ** 2 * sigmav
plt.plot(z_list, np.cumsum(n_ion[:, 0]+n_ion[:, 1])*Factor/(nb*0.76))
plt.xscale('log')
plt.yscale('log')
plt.show()