
#################################################################################################

# cosmo = Cosmology.setCosmology('planck1')
def CLfit(z):
    return 1.0+8e5*z**-1.4*np.exp(-(z/40)**2)

nb = 2.2e-7
# Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11
z_start = 1000.0
z_step = 0.01
z_list = np.array([0])
f = (Eg_list[1]/Eg_list[0])**1
double_shift = int(np.ceil(np.log(2.0)/np.log(f)))
for i in range(1000):
    z_list = np.append(z_list, f*(1.+z_list[-1])-1.0)
z_list = z_list[z_list <= z_start]
z_list = z_list[::-1]
#boost = np.interp(z_list, z_list_ann, boost_ann)
boost = CLfit(z_list)

import cosmolopy as cp
age_list = cp.distance.age(z_list,**cp.fidcosmo)

electrons = np.zeros(len(EgeV_list))
electrons = EgeV_list / 6.24e11
p2 = np.interp(EgeV_list, DM_spec_p[:,0]*1e9, DM_spec_p[:,1], left=0, right=0)
phot2 = np.interp(EgeV_list, DM_spec_phot[:,0]*1e9, DM_spec_phot[:,1], left=0, right=0)

# p2 = np.interp(EgeV_list, DM_spec_p[:,0]*1e9/4.0, DM_spec_p[:,1]/4.0, left=0, right=0)
# phot2 = np.interp(EgeV_list, DM_spec_phot[:,0]*1e9/4.0, DM_spec_phot[:,1]/4.0, left=0, right=0)


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
#Sigma_photoion = sigmaX(EgeV_list, 1, 1)*1e-18
Sigma_photoion[np.isnan(Sigma_photoion)] = 0.0


plt.plot(EgeV_list, Sigma_photoion)
plt.xscale('log'); plt.yscale('log')

n_ion = np.zeros([len(z_list), 10])
electronsN = np.zeros(len(EgeV_list))

for ii in range(len(z_list)-1):
    # plt.clf()
    # plt.subplot(211)
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
    N_naked = nb*c*tau*(1+z)**3
    N_H = N_naked*0.76
    N_He = N_naked*0.23/4
    N_O = N_naked*0.01/8
    # pair production
    temp_factor_H = (Sigma_photoion*N_H)
    temp_factor_H[EgeV_list < 13.6] = 0.0
    temp_factor_H[temp_factor_H>1.0] = 1.0
    temp = photons_particles * temp_factor_H
    photons_particles -= temp
    n_ion[ii, 0] += np.sum(temp*EgeV_list)/13.6/3

plt.plot(z_list, n_ion[:, 0])
plt.xscale('log')
plt.yscale('log')

sigmav = 2.16e-26
mx = 40e9

MatterDensityInGramsPerComovingCm = 5.7687745605587125e-30

Factor = (MatterDensityInGramsPerComovingCm / (mx * 1.783e-33)) ** 2 * sigmav
plt.plot(z_list, np.cumsum(n_ion[:, 0]+n_ion[:, 1])*Factor/(nb*0.76))
plt.plot(z_list, np.cumsum(n_ion[:, 4]+n_ion[:, 5])*Factor/(nb*0.23/4))
plt.plot(z_list, np.cumsum(n_ion[:, 8])*Factor/(nb*0.01/8))
plt.xscale('log')
# plt.yscale('log')

plt.plot(z_list, n_ion[:, 9])
plt.xscale('log')
plt.yscale('log')

plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2,'-')
plt.plot(EgeV_list, electronsN/np.gradient(Eg_list)*Eg_list**2,'--')
plt.xscale('log')
plt.yscale('log')
