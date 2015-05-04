import numpy as np
import matplotlib.pyplot as plt
from rstat.figurestyle import *
from rstat.crossections import sigmaH

# define_figure_style()
# http://arxiv.org/pdf/0807.1969.pdf
sigmaT = 6.6524e-25 #cm^2
me = 9.109e-28 # g
c = 2.9979e10 # cm/s
#hbar = 4.135668e-15/2.0/np.pi # eV*s
hbar = 1.0545e-27 # erg*s
kB = 1.380648e-16 # erg/K

def sigmakn(Eg, e, gamma):
    Gamma = 4*e*gamma/me/c**2
    eta = e*Eg/(me*c**2)**2
    q = Eg/Gamma/(gamma*me*c**2-Eg)
    # print np.sum((4*gamma**2)**-1 > q)
    # print np.sum(1 < q)
    G = 2.0*q*np.log(q)+(1.0+2.0*q)*(1.0-q)+2.0*eta*q*(1.0-q)
    G[G < 0] = 0
    G[(4*gamma**2)**-1 > q] = 0
    G[q > 1] = 0
    return 3.0*sigmaT/4.0/e/gamma**2*G


Eg_list = np.logspace(-2, 12, 1000)/6.24e11
temp = np.zeros([len(Eg_list), 4])
EgeV_list = Eg_list*6.24e11


from scipy.interpolate import InterpolatedUnivariateSpline


crossections = np.loadtxt('crossdata.dat')
crossectionsHe = np.loadtxt('crossdataHe.dat')
crossectionsO = np.loadtxt('crossdataO.dat')
crossectionsHe2 = np.loadtxt('photoionHe.dat')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-1]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '-', lw=1)
s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '-', lw=1)
plt.plot(EgeV_list, sigmaX(EgeV_list, 1, 1)*1e-18)
plt.plot(EgeV_list, sigmaX(EgeV_list, 8, 8)*1e-18)
plt.xscale('log')
plt.yscale('log')

plt.plot(EgeV_list, sigmaX(EgeV_list, 1, 1) / sigmaX(EgeV_list, 8, 8))
plt.xscale('log')
plt.yscale('log')


plt.plot(EgeV_list/1e9, 6.3e-18 * (EgeV_list/13.6)**-3, '-', lw=1)
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
factor = 10**np.interp(np.log10(EgeV_list), np.log10(crossectionsHe2[:,0]), np.log10(crossectionsHe2[:,1]))
plt.plot(EgeV_list/1e9, factor*1e-18, '-', lw=1)
plt.xscale('log')
plt.yscale('log')

factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,-1]))
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,1]))
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,2]))
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,3]))
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,5]))
plt.plot(EgeV_list/1e9, factor*1e-24, '-', lw=1)
plt.xscale('log')
plt.yscale('log')


# def nph(e, T):
#     ec = kB*T
#     return 1.0/(np.pi**2*hbar**3*c**3)*e**2/(np.exp(e/ec)-1.0)
#
# def Int1(Eg, gamma, T):
#     e_list = np.logspace(np.log10(kB*T)-4, np.log10(kB*T)+3, 100)
#     # e_list = np.logspace(np.log10(kB*T)-0.0001, np.log10(kB*T)+0.0001, 10)
#     nph_list = nph(e_list, T)
#     sigmakn_list = sigmakn(Eg, e_list, gamma)
#     return np.trapz(nph_list*sigmakn_list, e_list)
#
# Int1vec = np.vectorize(Int1)
#
# def dNdg(gamma):
#     E = gamma*me*c**2
#     EeV = E*6.24e11
#     E2dNdE = 10**np.interp(np.log10(EeV), np.log10(px[:, 0]*1e9), np.log10(2*px[:, 1]), left=-100, right=-100)
#     return E2dNdE/E**2/(me*c**2)
#
# def dNdgD(gamma):
#     E = gamma*me*c**2
#     EeV = E*6.24e11
#     E2dNdE = np.exp(-(gamma-6653)**2/1e3)*1e20
#     # E2dNdE = np.array([0,1,0])*1e20
#     return E2dNdE/E**2/(me*c**2)
#
# def P(Eg, T):
#     gamma_list = np.logspace(-2, 6, 300)
#     # gamma_list = np.array([100,6653,10000])
#     Ngamma_list = dNdgD(gamma_list)
#     Int1_list = Int1vec(Eg, gamma_list, T)
#     # Int1_list /= np.sum(Int1_list)
#     temp = Ngamma_list*Int1_list
#     # print np.sum(Int1_list<0), np.sum(Int1_list>0)
#     temp[np.isnan(temp)] = 0
#     result = c / Eg * np.trapz(temp * gamma_list, gamma_list)
#     return result
#
# Pvec = np.vectorize(P)


Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11

z = 0.0
T_CMB = 2.7*(1.0+z)
nb = (1.0+z)**3*2.2e-7
# T_list=np.array([1000])

photons = np.zeros(len(Eg_list))
# photons[np.where(np.abs(Eg_list-T_CMB/11500) == np.abs(Eg_list-T_CMB/11500).min())[0][0]] = 511
CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)
electrons = np.zeros(1)
electrons[0] = 1e8 / 6.24e11
# electrons /= np.trapz(electrons, Eg_list)
# print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)

tau = 3.14e7*1e4
photons_particles = np.zeros(len(Eg_list))

# Do electrons-photons interactions
for j in range(100000):
    chances = np.random.random([len(electrons), 2])
    for i in range(len(electrons)):
        gamma = electrons[i] / (me*c**2)
        sigma_IC = np.zeros(len(Eg_list))
        for j1 in range(len(photons)):
            sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
            sigma_IC[j1] = np.trapz(sigma_temp, Eg_list)
        sigma_IC_total = sigma_IC*CMBphotons*np.gradient(Eg_list)
        probabilities = sigma_IC_total * c * tau
        probabilities = np.cumsum(probabilities)
        j1 = np.where(chances[i, 0] > probabilities)[0]
        if len(j1)<len(probabilities):
            j1 = j1[-1]
            sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
            sigma_temp = np.cumsum(sigma_temp)
            sigma_temp /= sigma_temp.max()
            # print j2
            j2 = np.where(chances[i, 1] > sigma_temp)[0][-1]
            photons_particles[j2] += 1
            electrons[i] -= Eg_list[j2]
    if j%100==0:
        print electrons[0]


plt.plot(EgeV_list, photons_particles*Eg_list**2)
plt.xscale('log')
plt.yscale('log')

####################################

z=100
T_CMB=2.7*(1+z)
CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
plt.plot(EgeV_list, CMBphotons*Eg_list**2)
z=10
T_CMB=2.7*(1+z)
CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
plt.plot(EgeV_list, CMBphotons*Eg_list**2)
plt.xscale('log')
plt.yscale('log')


####################################

electrons = np.zeros(1)
electrons[0] = 1e8 / 6.24e11
# electrons /= np.trapz(electrons, Eg_list)
# print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)
z=100
T_CMB=2.7*(1+z)

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

tau = 3.14e7*1e2
photons_particles = np.zeros(len(Eg_list))

# Do electrons-photons interactions
for j in range(100):
    chances = np.random.random([len(electrons), 2])
    for i in range(len(electrons)):
        gamma = electrons[i] / (me*c**2)
        sigma_IC = np.zeros(len(Eg_list))
        for j1 in range(len(photons)):
            sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
            # sigma_IC[j1] = np.trapz(sigma_temp, Eg_list)
            # sigma_IC_total = sigma_IC[j1]*CMBphotons[j1]*np.gradient(Eg_list)[j1]
            # probability = sigma_IC_total * c * tau
            temp = sigma_temp * c * tau * CMBphotons[j1] * np.gradient(Eg_list)[j1]
            photons_particles += temp
            # print np.sum(temp),
            probability = np.sum(temp*Eg_list)
            electrons[i] -= probability
            # print probability
            # sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
            # sigma_temp = np.cumsum(sigma_temp)
            # sigma_temp /= sigma_temp.max()
            # # print j2
            # j2 = np.where(chances[i, 1] > sigma_temp)[0][-1]
            # photons_particles[j2] += 1
            # electrons[i] -= Eg_list[j2]
        print np.log10(electrons[0]*6.24e11)
    if j % 200 == 199:
        plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)


print np.log10(np.sum(photons_particles*Eg_list*6.24e11))
print np.log10(np.sum(photons_particles*Eg_list*6.24e11)+electrons[0]*6.24e11)
plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)
plt.plot(EgeV_list, CMBphotons*Eg_list**2)
plt.xscale('log')
plt.yscale('log')

plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)
plt.plot(EgeV_list, (EgeV_list/1e9)**0.5*7e-2 *(np.exp(-(EgeV_list/1e9*6e3*100))))
plt.xscale('log')
plt.yscale('log')

def ICon30CMB(EgeV_list, E0, z):
    return (EgeV_list/1e9)**0.5*5e-2 *(np.exp(-(EgeV_list/1e9*3.4e3*(1e9/E0/((1.+z)/61)**0.5)**2)))

plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)
plt.plot(EgeV_list, ICon30CMB(EgeV_list, 1e8, z))
plt.xscale('log')
plt.yscale('log')

plt.plot(p[:,0], p[:,1]/p[:,0]**2*np.gradient(p[:,0]),'--',label=r'$\mathrm{positrons}$')
plt.xscale('log')
plt.yscale('log')

spec = np.zeros(len(EgeV_list))
for i in range(len(p[:,0])):
    N = (p[:,1]/p[:,0]**2*np.gradient(p[:,0]))[i]
    spec += N*ICon30CMB(EgeV_list, p[i,0]*1e9)


# plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)

factor = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,-1]))
plt.plot(EgeV_list, 1/3e10/(factor*1e-24*nb*(1+100)**3)/3e7, '-', lw=1)
plt.plot(EgeV_list, spec)
plt.xscale('log')
plt.yscale('log')


#################################################################################################
#
# # Eg_list = np.logspace(-6, 12, 1000)/6.24e11
# nu_list = Eg_list / (hbar*2*np.pi)
# EgeV_list = Eg_list*6.24e11
# z_start = 1000.0
# z_step = 0.01
# z_list = np.array([0])
# f = (Eg_list[1]/Eg_list[0])**1
# double_shift = int(np.ceil(np.log(2.0)/np.log(f)))
# for i in range(1000):
#     z_list = np.append(z_list, f*(1.+z_list[-1])-1.0)
# z_list = z_list[z_list <= z_start]
# z_list = z_list[::-1]
#
# import cosmolopy as cp
# age_list = cp.distance.age(z_list,**cp.fidcosmo)
#
# electrons = np.zeros(len(EgeV_list))
# electrons = EgeV_list / 6.24e11
# p2 = np.interp(EgeV_list, p[:,0]*1e9, p[:,1], left=0, right=0)
# phot2 = np.interp(EgeV_list, phot[:,0]*1e9, phot[:,1], left=0, right=0)
#
# # plt.plot(p[:,0], p[:,1]); plt.plot(EgeV_list/1e9, p2);
# # plt.xscale('log');
# # plt.yscale('log')
#
# electronsN = (p2/(EgeV_list/1e9)**2*np.gradient(EgeV_list/1e9))
#
# # electrons = np.array([1e8/6.24e11])
# # electrons /= np.trapz(electrons, Eg_list)
# # print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)
# z = z_list[0]
# T_CMB = 2.7*(1+z)
# CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
# print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)
#
# tau = 3.14e7*1e2
# photons_particles = np.zeros(len(Eg_list))
# photons_particles += (phot2/(EgeV_list/1e9)**2*np.gradient(EgeV_list/1e9))
# temp = 0
#
# Sigma_pair_prod = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,4]+crossections[:,5]))*1e-24
# Sigma_pair_prod[np.isnan(Sigma_pair_prod)] = 0.0
#
# Sigma_photoion = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,3]))*1e-24
# Sigma_photoion[np.isnan(Sigma_photoion)] = 0.0
#
# Sigma_photoion_He = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,3]))*1e-24
# Sigma_photoion_He[np.isnan(Sigma_photoion_He)] = 0.0
#
# Sigma_photoion_O = sigmaX(EgeV_list, 8, 8)*1e-18
#
# Sigma_collisional_ion_He = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,2]))*1e-24
# Sigma_collisional_ion_He[np.isnan(Sigma_collisional_ion_He)] = 0.0
#
# Sigma_collisional_scatter = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,1]))*1e-24
# Sigma_collisional_scatter[np.isnan(Sigma_collisional_scatter)] = 0.0
#
# n_ion = np.zeros([len(z_list), 10])
#
# for ii in range(len(z_list)-1):
#     plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2,'-')
#     plt.plot(EgeV_list, electronsN/np.gradient(Eg_list)*Eg_list**2,'--')
#     plt.xscale('log')
#     plt.yscale('log')
#     print ii
#     dt = age_list[ii+1]-age_list[ii]
#     tau = dt
#     print tau
#     z = 0.5 * (z_list[ii+1] + z_list[ii])
#     T_CMB = 2.7*(1+z)
#     # redshift the system
#     photons_particles[:-1] = photons_particles[1:]
#     photons_particles[-1] = 0.
#     # update CMB
#     CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
#     # add IC photons
#     for i in range(len(electrons)):
#         photons_particles += electronsN[i] * ICon30CMB(EgeV_list, electrons[i]*6.24e11, z) * np.gradient(Eg_list) / Eg_list**2
#     electronsN *= 0
#     # column density
#     N_naked = nb*c*tau*(1+z)**3
#     N_H = N_naked*0.76
#     N_He = N_naked*0.23/4
#     N_O = N_naked*0.01/8
#     # pair production
#     temp_factor = (Sigma_pair_prod*N_naked)
#     temp_factor[temp_factor > 1.0] = 1.0
#     temp = photons_particles * temp_factor
#     photons_particles -= temp
#     electronsN[:-double_shift] += temp[double_shift:]*2
#     # photoion
#     temp_factor_H = (Sigma_photoion*N_H)
#     temp_factor_He = (Sigma_photoion_He*N_He)
#     temp_factor_O = (Sigma_photoion_O*N_O)
#     temp_factor = temp_factor_H + temp_factor_He + temp_factor_O
#     temp_factor_H[EgeV_list < 13.6] = 0.0
#     temp_factor_He[EgeV_list < 24.6] = 0.0
#     temp_factor_H[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
#     temp_factor_He[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
#     temp_factor_O[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
#     temp_factor[temp_factor > 1.0] = 1.0
#     temp = photons_particles * temp_factor
#     photons_particles -= temp
#     n_ion[ii, 0] += np.sum(photons_particles * temp_factor_H)
#     n_ion[ii, 4] += np.sum(photons_particles * temp_factor_He)
#     n_ion[ii, 8] += np.sum(photons_particles * temp_factor_O)
#     # collisional ion
#     temp_factor = (Sigma_collisional_ion*N_naked)
#     temp_factor[temp_factor > 1.0] = 1.0
#     temp_factor[EgeV_list < 13.6] = 0.0
#     temp = photons_particles * temp_factor
#     photons_particles -= temp
#     n_ion[ii, 1] += np.sum(temp)
#
# plt.plot(z_list, n_ion[:, 0])
# plt.plot(z_list, n_ion[:, 4])
# plt.xscale('log')
# plt.yscale('log')
#
# plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2,'-')
# plt.plot(EgeV_list, electronsN/np.gradient(Eg_list)*Eg_list**2,'--')
# plt.xscale('log')
# plt.yscale('log')
#

#################################################################################################

cosmo = Cosmology.setCosmology('planck1')

z_list_ann = z_list.copy()
boost_ann = 1+rho_dm_2[:]

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
boost = np.interp(z_list, z_list_ann, boost_ann)

import cosmolopy as cp
age_list = cp.distance.age(z_list,**cp.fidcosmo)

electrons = np.zeros(len(EgeV_list))
electrons = EgeV_list / 6.24e11
p2 = np.interp(EgeV_list, p[:,0]*1e9, p[:,1], left=0, right=0)
phot2 = np.interp(EgeV_list, phot[:,0]*1e9, phot[:,1], left=0, right=0)

p2 = np.interp(EgeV_list, p[:,0]*1e9/4.0, p[:,1]/4.0, left=0, right=0)
phot2 = np.interp(EgeV_list, phot[:,0]*1e9/4.0, phot[:,1]/4.0, left=0, right=0)


z = z_list[0]
T_CMB = 2.7*(1+z)
CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

tau = 3.14e7*1e2
photons_particles = np.zeros(len(Eg_list))
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

Sigma_photoion_O = sigmaX(EgeV_list, 8, 8)*1e-18

Sigma_collisional_ion_He = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,2]))*1e-24
Sigma_collisional_ion_He[np.isnan(Sigma_collisional_ion_He)] = 0.0

Sigma_collisional_scatter = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,1]))*1e-24
Sigma_collisional_scatter[np.isnan(Sigma_collisional_scatter)] = 0.0

# s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
# Sigma_photoion_O = 10**s(np.log10(EgeV_list/1e6))*1e-24
# Sigma_photoion_O[np.isnan(Sigma_photoion_He)] = 0.0

n_ion = np.zeros([len(z_list), 10])
electronsN = np.zeros(len(EgeV_list))

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
    photons_particles += (phot2/(EgeV_list/1e9)**2*np.gradient(EgeV_list/1e9))*boost[ii]*(1.0+z)**3*tau
    electronsN = 2*(p2/(EgeV_list/1e9)**2*np.gradient(EgeV_list/1e9))*boost[ii]*(1.0+z)**3*tau
    n_ion[ii,9] = (1.0+z)**3*tau*boost[ii]
    # add IC photons
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
    temp_factor_O = (Sigma_photoion_O*N_O)
    temp_factor_H[EgeV_list < 13.6] = 0.0
    temp_factor_He[EgeV_list < 24.6] = 0.0
    temp_factor_O[EgeV_list < 13.6] = 0.0
    temp_factor = temp_factor_H + temp_factor_He + temp_factor_O
    temp_factor_H[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor_He[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor_O[temp_factor > 1.0] /= temp_factor[temp_factor > 1.0]
    temp_factor[temp_factor > 1.0] = 1.0
    temp = photons_particles * temp_factor
    photons_particles -= temp
    n_ion[ii, 0] += np.sum(photons_particles * temp_factor_H)
    n_ion[ii, 4] += np.sum(photons_particles * temp_factor_He)
    n_ion[ii, 8] += np.sum(photons_particles * temp_factor_O)
    # plt.subplot(212)
    # plt.plot(EgeV_list, temp_factor_H, '-', lw=2)
    # plt.plot(EgeV_list, temp_factor_He, '-', lw=1)
    # plt.xscale('log')
    # plt.yscale('log')
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

plt.plot(z_list, n_ion[:, 0])
plt.plot(z_list, n_ion[:, 4])
plt.plot(z_list, n_ion[:, 8])
plt.xscale('log')
plt.yscale('log')

sigmav = 2.16e-26
mx = 4e9
nb = 2.2e-7

Factor = ((cosmo.matterDensity(0.0) / Cosmology.AST_kpc**3 * Cosmology.AST_Msun) / (mx * 1.783e-33)) ** 2 * sigmav
plt.plot(z_list, np.cumsum(n_ion[:, 0]+n_ion[:, 1])*Factor/(nb*0.76))
# plt.xscale('log')
# plt.yscale('log')

Factor = ((cosmo.matterDensity(0.0) / Cosmology.AST_kpc**3 * Cosmology.AST_Msun) / (mx * 1.783e-33)) ** 2 * sigmav
plt.plot(z_list, np.cumsum(n_ion[:, 4]+n_ion[:, 5])*Factor/(nb*0.23/4))
# plt.xscale('log')
# plt.yscale('log')

Factor = ((cosmo.matterDensity(0.0) / Cosmology.AST_kpc**3 * Cosmology.AST_Msun) / (mx * 1.783e-33)) ** 2 * sigmav
plt.plot(z_list, np.cumsum(n_ion[:, 8])*Factor/(nb*0.01/8))
plt.xscale('log')
plt.yscale('log')

plt.plot(z_list, n_ion[:, 9])
plt.xscale('log')
plt.yscale('log')

plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2,'-')
plt.plot(EgeV_list, electronsN/np.gradient(Eg_list)*Eg_list**2,'--')
plt.xscale('log')
plt.yscale('log')

