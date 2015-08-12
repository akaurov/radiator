import numpy as np
import cosmolopy
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

# Do electrons-photons interactions
E0_list = np.logspace(1, 12, 250)
z = 50.0
Eth_IC = 10.
delta = 1000000.0
T_CMB = 2.73*(1+z)
xiH = 0.1
nb = cosmolopy.cden.baryon_densities(**cosmolopy.fidcosmo)[0] / cosmolopy.constants.Mpc_cm**3 * cosmolopy.constants.M_sun_g / cosmolopy.constants.m_p_g * 0.04
nb *= (1.0+delta)
T = 10**4
precision = 0.01
accurate_cmb = True

mass_abundance = np.zeros([10,10])
mass_abundance[ 1, 1] = 0.754 * (1.0 - xiH)
mass_abundance[ 1, 0] = 0.754 * xiH
mass_abundance[ 2, 2] = 0.246 * (1.0 - xiH)
mass_abundance[ 2, 1] = 0.246 * xiH

nbHI = mass_abundance[ 1, 1] * nb
nbHeI = mass_abundance[ 2, 2] * nb / 4.0
nbHeII = mass_abundance[ 2, 1] * nb / 4.0
nbe = mass_abundance[ 1, 0] * nb + mass_abundance[ 2, 1] * nb / 4.0 + mass_abundance[ 2, 0] * nb / 4.0

Eg_list = np.logspace(-6, np.log10(E0_list.max()), 300)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
N_CMB = np.trapz(CMBphotons, Eg_list)
E_CMB_av = np.trapz(CMBphotons*Eg_list, Eg_list) / np.trapz(CMBphotons, Eg_list)
print 'number of photons per cm^3: ', N_CMB
print 'average CMB photon energy [eV]: ', E_CMB_av *6.24e11

# tau_list = np.zeros([len(E0_list), 10000]) * np.nan
# IC_total_list = np.zeros([len(E0_list), 10000])
# coll_ion_list = np.zeros([len(E0_list), 10000])
# coll_ex_list = np.zeros([len(E0_list), 10000])
# coll_loss_list = np.zeros([len(E0_list), 10000])

tau = 3.14e7*1e-5
T_tau = 0
photons_particles_all = np.zeros([len(Eg_list), len(E0_list)])
coll_ion_N = 0
coll_ion = 0
coll_ion = 0
coll_ion = 0
coll_loss = 0
coll_ex = 0
IC_total = 0
# E0_list = np.array([1e6])
results = np.zeros([len(E0_list), 10])
MC_N_list = np.zeros([len(E0_list)])

for i_E in range(len(E0_list)):
    E0 = E0_list[i_E]

    if E0 < 1e2:
        MC_N = 1000
        precision = 0.01
    elif E0 < 1e4:
        MC_N = 100
        precision = 0.01
    else:
        MC_N = 10
        precision = 0.01

    MC_N_list[i_E] = MC_N

    coll_ion = 0
    coll_ion_HeI = 0
    coll_ion_HeII = 0
    coll_loss = 0
    coll_ex = 0
    coll_ex_HeI = 0
    coll_ex_HeII = 0
    IC_total = 0
    IC_soft = 0
    IC_hard = 0
    photons_particles = np.zeros(len(Eg_list))

    gamma = (E0 / 6.24e11 + (me*c**2)) / (me*c**2)
    v = c*np.sqrt(1.-1./gamma**2)
    sigma_IC = np.zeros(len(Eg_list))
    total_probability = 0
    photons_particles_total_add = np.zeros(len(Eg_list))
    tau_IC = 3.14e100
    if E0 > Eth_IC:
        # accurate_cmb = electrons[i] < 1/6.24e11
        if accurate_cmb:
            for j1 in range(len(Eg_list)):
                sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                temp = sigma_temp * v * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                probability = np.sum(temp*Eg_list)
                total_probability += probability
                photons_particles_total_add += temp
        else:
            sigma_temp = sigmakn(Eg_list, E_CMB_av, gamma) * np.gradient(Eg_list)
            temp = sigma_temp * v * N_CMB
            probability = np.sum(temp*Eg_list)
            total_probability += probability
            photons_particles_total_add += temp
        tau_IC = E0 / 6.24e11 / total_probability

    for iii in range(MC_N):
        electrons = np.zeros(1)
        electrons[0] = E0 / 6.24e11
        j = 0
        while (electrons.sum()*6.24e11 > 1e-10):
            j = j+1
            # chances = np.random.random([len(electrons), 2])
            for i in range(len(electrons)):
                if electrons[i] > 10.2/6.24e11:
                    gamma = (electrons[i] + (me*c**2)) / (me*c**2)
                    v = c*np.sqrt(1.-1./gamma**2)
                    sigma_IC = np.zeros(len(Eg_list))

                    eedEdt_now = eedEdt(electrons[i], nbe*(1+z)**3, T)
                    tau_ex = nbHI*(1+z)**3*sigmaHex(electrons[i])*v
                    tau_ion = nbHI*(1+z)**3*sigmaHion(electrons[i])*v
                    tau_ex_HeI = nbHeI*(1+z)**3*sigmaHe(electrons[i], 'ex', 'I')*v
                    tau_ion_HeI = nbHeI*(1+z)**3*sigmaHe(electrons[i], 'ion', 'I')*v
                    tau_ex_HeII = nbHeII*(1+z)**3*sigmaHe(electrons[i], 'ex', 'II')*v
                    tau_ion_HeII = nbHeII*(1+z)**3*sigmaHe(electrons[i], 'ion', 'II')*v

                    if electrons[i] > 1e5/6.24e11:
                        # Approximate mode
                        # tau = precision*min([1./tau_ex, 1./tau_ion, tau_IC, electrons[i]/eedEdt_now])
                        secondary_energy = rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                        energy_to_distribute = electrons[i] - E0_list[i_E-1]/6.24e11 + 1/6.24e11
                        e_to_ex = tau_ex * 10.2 / 6.24e11
                        e_to_ion = tau_ion * 13.6/ 6.24e11
                        e_to_sec = tau_ion * secondary_energy
                        e_to_ee = eedEdt_now
                        e_to_ic = total_probability
                        e_factor = energy_to_distribute / (e_to_ee+e_to_ex+e_to_ion+e_to_ic+e_to_sec)
                        e_to_ex *= e_factor
                        e_to_ion *= e_factor
                        e_to_ee *= e_factor
                        e_to_ic *= e_factor
                        e_to_sec *= e_factor
                        coll_ex += e_to_ex
                        coll_loss += e_to_ee
                        IC_total += e_to_ic
                        coll_ion += e_to_ion
                        photons_particles += photons_particles_total_add*e_factor
                        e_factor_int = np.floor(e_to_sec / secondary_energy)
                        print electrons[i]*6.24e11, energy_to_distribute
                        electrons[i] -= energy_to_distribute
                        electrons = np.append(electrons, np.ones(e_factor_int)*e_to_sec/e_factor_int)
                    else:
                        if electrons[i] > 1e5/6.24e11:
                            tau = precision*min([1./tau_ex, 1./tau_ion, tau_IC, electrons[i]/eedEdt_now])
                        elif electrons[i] > 13.6/6.24e11:
                            tau = precision*min([1./tau_ex, 1./tau_ion, tau_IC, electrons[i]/eedEdt_now])
                        else:
                            tau = precision*min([1./tau_ex, tau_IC, electrons[i]/eedEdt_now])
                        T_tau += tau
                        e0z = electrons[i].copy()
                        electrons[i] -= total_probability*tau
                        IC_total += total_probability*tau
                        photons_particles += photons_particles_total_add*tau
                        # if electrons[i] > 1e5/6.24e11:
                        #     print total_probability, np.sum(photons_particles_total_add*Eg_list)
                        rand_a = np.random.rand(10)
                        if electrons[i] > 10.2 / 6.24e11:
                            if rand_a[0] < (1.0 - np.exp(-tau_ex*tau)):
                                electrons[i] -= 10.2 / 6.24e11
                                coll_ex += 10.2 / 6.24e11
                        if electrons[i] > 20.0 / 6.24e11:
                            if rand_a[1] < (1.0 - np.exp(-tau_ex_HeI*tau)):
                                electrons[i] -= 20.0 / 6.24e11
                                coll_ex += 20.0 / 6.24e11
                        if electrons[i] > 41.0 / 6.24e11:
                            if rand_a[2] < (1.0 - np.exp(-tau_ex_HeII*tau)):
                                electrons[i] -= 41.0 / 6.24e11
                                coll_ex += 41.0 / 6.24e11
                        if electrons[i] > 13.6 / 6.24e11:
                            if rand_a[3] < (1.0 - np.exp(-tau_ex*tau)):
                                temp = rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= 13.6 / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion += 13.6 / 6.24e11
                        if electrons[i] > 24.6 / 6.24e11:
                            if rand_a[4] < (1.0 - np.exp(-tau_ex*tau)):
                                temp = rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= 24.6 / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion += 24.6 / 6.24e11
                        if electrons[i] > 54.4 / 6.24e11:
                            if rand_a[5] < (1.0 - np.exp(-tau_ex*tau)):
                                temp = rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= 54.4 / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion += 54.4 / 6.24e11
                        if electrons[i] > eedEdt_now * tau:
                            electrons[i] -= eedEdt_now * tau
                            coll_loss += eedEdt_now * tau
                        # tau_list[iii, j] = T_tau/3.14e7
                        # print i, j, electrons[i]/E0*6.24e11
                        # print IC_total*6.24e11, coll_ion*6.24e11, coll_ex*6.24e11, coll_loss*6.24e11
                else:
                    coll_loss += electrons[i]
                    electrons[i] = 0
            # interpolation mode
            # inter_mask = np.where(electrons < E0_list[i_E-1]/6.24e11)[0]
            inter_mask = np.where((electrons < E0_list[i_E-1]/6.24e11))[0]# & (electrons < Eth_IC/6.24e11))[0]
            IC_total_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 0])
            IC_soft_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 4])
            IC_hard_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 5])
            coll_ion_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 1])
            coll_ex_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 2])
            coll_loss_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 3])
            IC_total += np.sum(IC_total_frac * electrons[inter_mask])
            IC_soft += np.sum(IC_soft_frac * electrons[inter_mask])
            IC_hard += np.sum(IC_hard_frac * electrons[inter_mask])
            coll_ion += np.sum(coll_ion_frac * electrons[inter_mask])
            coll_ex += np.sum(coll_ex_frac * electrons[inter_mask])
            coll_loss += np.sum(coll_loss_frac * electrons[inter_mask])
            electrons = np.delete(electrons, inter_mask)
    results[i_E,:6] = IC_total*6.24e11, coll_ion*6.24e11, coll_ex*6.24e11, coll_loss*6.24e11, \
                      (IC_soft+np.sum((photons_particles*Eg_list)[EgeV_list < 13.6]))*6.24e11, \
                      (IC_hard+np.sum((photons_particles*Eg_list)[EgeV_list >= 13.6]))*6.24e11

    photons_particles_all[:, i_E] = photons_particles
    results[i_E, :] /= MC_N*E0_list[i_E]
    print E0, results[i_E,:6], np.sum(results[i_E,:4])


np.savez('output/%05i-%01.5f-%08.1f-snap.npz'%(z,xiH,delta), E0_list=E0_list, results=results, photons_particles_all=photons_particles_all, Eg_list=Eg_list, EgeV_list=EgeV_list)

plt.plot(E0_list, np.sum(results[:, 1:4], 1))
plt.plot(E0_list, results[:, 4])
plt.plot(E0_list, results[:, 5])
# plt.plot(E0_list, results[:, 1])
# plt.plot(E0_list, results[:, 2])
# plt.plot(E0_list, results[:, 3])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$E\;\mathrm{[eV]}$')


    for i in range(len(E_ph_interaction_list)):
        results[i_E, 7] +=

# plt.savefig('%05i-%1.5f-%5.1f-fig.png'%(z,xi,delta), dpi=300)
# # plt.yscale('log')
#
#
# plt.plot(EgeV_list, 3*(13.6+8)*(1-np.exp(-nb*(1+z)**3*sigmaHex(Eg_list)*c*np.sqrt(1-511e3**2/(511e3+EgeV_list)**2))))
# plt.plot(EgeV_list, 2.3e-13*1e9*(1+z)**4*(EgeV_list/100e9)**2)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# j += 1
# IC_total_list[iii, j] = 0
# coll_ion_list[iii, j] = 0
# coll_ex_list[iii, j] = 0
# coll_loss_list[iii, j] = electrons[i]
#
#         # if j % 200 == 199:
#             # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)
#
# plt.clf()
# total_energy = IC_total_list+coll_ion_list+coll_ex_list+coll_loss_list
# plt.plot(E0_list, IC_total_list.sum(1)/total_energy.sum(1), '-k')
# plt.plot(E0_list, coll_ion_list.sum(1)/total_energy.sum(1), '--k')
# plt.plot(E0_list, coll_ex_list.sum(1)/total_energy.sum(1), '-.k')
# plt.plot(E0_list, coll_loss_list.sum(1)/total_energy.sum(1), ':k')
# plt.xscale('log')
# plt.xlim([10,1e9])
# # plt.yscale('log')
# plt.xlabel(r'$E_0\;\mathrm{[eV]}$')
# plt.ylabel(r'$\mathrm{Fractional\;contribution}$')
# fixlogax(plt.gca(), 'x')
# # plt.savefig('figures/FracVsE0z10xi001.pdf')

#########################