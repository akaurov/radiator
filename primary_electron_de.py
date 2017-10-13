#!python

import sys

import numpy as np
import cosmolopy
import matplotlib.pyplot as plt

from multiprocessing import Pool
n_proc = 4

from radiator.figurestyle import *
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
# E0_list = np.logspace(1, 12, 250)
E0_list = np.logspace(1, 12, 1000)
E_extr_thresh = 1e12
z = float(sys.argv[1])
Eth_IC = 1e3
delta = float(sys.argv[2])
T_CMB = 2.73*(1+z)
xiH = float(sys.argv[3])
nb = cosmolopy.cden.baryon_densities(**cosmolopy.fidcosmo)[0] / cosmolopy.constants.Mpc_cm**3 * cosmolopy.constants.M_sun_g / cosmolopy.constants.m_p_g * 0.04
nb *= (1.0+delta)
T = 4e4
precision = 0.01
accurate_cmb = False
cs = cross_sections(cs={'photion': sys.argv[4], 'collion': sys.argv[5], 'collex': sys.argv[6]})

mass_abundance = np.zeros([10,10])
mass_abundance[ 1, 1] = 0.754 * (1.0 - xiH)
mass_abundance[ 1, 0] = 0.754 * xiH
mass_abundance[ 2, 2] = 0.246 * (1.0 - xiH)
mass_abundance[ 2, 1] = 0.246 * xiH

nbHI = mass_abundance[ 1, 1] * nb
nbHeI = mass_abundance[ 2, 2] * nb / 4.0
nbHeII = mass_abundance[ 2, 1] * nb / 4.0
nbe = mass_abundance[ 1, 0] * nb + mass_abundance[ 2, 1] * nb / 4.0 + mass_abundance[ 2, 0] * nb / 4.0

Eg_list = np.logspace(-6, np.log10(E0_list.max()), 100)/6.24e11
EgeV_list = Eg_list*6.24e11

nu_list = Eg_list / (hbar*2*np.pi)

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
N_CMB = np.trapz(CMBphotons, Eg_list)
E_CMB_av = np.trapz(CMBphotons*Eg_list, Eg_list) / np.trapz(CMBphotons, Eg_list)
print 'number of photons per cm^3: ', N_CMB
print 'average CMB photon energy [eV]: ', E_CMB_av * 6.24e11

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

results = np.zeros(len(E0_list), dtype=([('H_I_ion', 'f4'), ('H_I_ex', 'f4'),
                                         ('He_I_ion', 'f4'), ('He_I_ex', 'f4'),
                                         ('He_II_ion', 'f4'), ('He_II_ex', 'f4'),
                                         ('ee', 'f4'),
                                         ('IC', 'f4')]))
MC_N_list = np.zeros([len(E0_list)])



for i_E in range(0, len(E0_list)):
    print i_E
    E0 = E0_list[i_E]

    if E0 < 1e2:
        MC_N = 1000
        precision = 0.01
    elif E0 < 1e5:
        MC_N = 100
        precision = 0.01
    else:
        MC_N = 1
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
    tau_IC = 3.14e100
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

                    eedEdt_now = cs.eedEdt(electrons[i], nbe*(1+z)**3, T)
                    # tau_ex = nbHI*(1+z)**3*sigmaHex(electrons[i])*v
                    # tau_ion = nbHI*(1+z)**3*sigmaHion(electrons[i])*v
                    tau_ion = nbHI*(1+z)**3*cs.HI_ion_e(electrons[i]*6.24e11)*v
                    tau_ex = nbHI*(1+z)**3*cs.HI_ex_e(electrons[i]*6.24e11)*v
                    # tau_ex = nbHI*(1+z)**3*sigmaRBEQ(electrons[i], 13.6057, 13.6057, 1)*v/0.85
                    # tau_ion = nbHI*(1+z)**3*sigmaRBEQ(electrons[i], 13.6, 13.6, 1)*v
                    tau_ex_HeI = nbHeI*(1+z)**3*cs.HeI_ion_e(electrons[i]*6.24e11)*v
                    tau_ion_HeI = nbHeI*(1+z)**3*cs.HeI_ex_e(electrons[i]*6.24e11)*v
                    tau_ex_HeII = nbHeII*(1+z)**3*cs.HeII_ion_e(electrons[i]*6.24e11)*v
                    tau_ion_HeII = nbHeII*(1+z)**3*cs.HeII_ex_e(electrons[i]*6.24e11)*v

                    sigma_IC = np.zeros(len(Eg_list))
                    total_probability = 0
                    photons_particles_total_add = np.zeros(len(Eg_list))

                    if electrons[i]*6.24e11 > Eth_IC:
                        # accurate_cmb = electrons[i] < 1/6.24e11
                        if accurate_cmb:
                            for j1 in range(len(Eg_list)):
                                sigma_temp = cs.sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                                temp = sigma_temp * v * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                                probability = np.sum(temp*Eg_list)
                                total_probability += probability
                                photons_particles_total_add += temp
                        else:
                            sigma_temp = cs.sigmakn(Eg_list, E_CMB_av, gamma) * np.gradient(Eg_list)
                            temp = sigma_temp * v * N_CMB
                            probability = np.sum(temp*Eg_list)
                            total_probability += probability
                            photons_particles_total_add += temp

                    if electrons[i] > E_extr_thresh/6.24e11:
                        # Approximate mode
                        # tau = precision*min([1./tau_ex, 1./tau_ion, tau_IC, electrons[i]/eedEdt_now])
                        rhoE_list = cs.rhoE(np.ones(1000)*electrons[i], 8.0/6.24e11)
                        secondary_energy = np.mean(rhoE_list)
                        energy_to_distribute = electrons[i] * 0.3
                        e_to_ion = tau_ion * cs.const_HI_ion / 6.24e11
                        e_to_ex = tau_ex * cs.const_HI_ex / 6.24e11
                        e_to_ion_HeI = tau_ion_HeI * cs.const_HeI_ion / 6.24e11
                        e_to_ex_HeI = tau_ex_HeI * cs.const_HeI_ex / 6.24e11
                        e_to_ion_HeII = tau_ion_HeII * cs.const_HeII_ion / 6.24e11
                        e_to_ex_HeII = tau_ex_HeII * cs.const_HeII_ex / 6.24e11
                        e_to_sec = (tau_ion+tau_ion_HeI+tau_ion_HeII) * secondary_energy
                        e_to_ee = eedEdt_now
                        e_to_ic = total_probability
                        e_factor = energy_to_distribute / (e_to_ee+e_to_ex+e_to_ion+e_to_ic+e_to_sec+
                                                           e_to_ion_HeI+e_to_ion_HeII+
                                                           e_to_ex_HeI+e_to_ex_HeII)
                        e_to_ex *= e_factor
                        e_to_ion *= e_factor
                        e_to_ex_HeI *= e_factor
                        e_to_ion_HeI *= e_factor
                        e_to_ex_HeII *= e_factor
                        e_to_ion_HeII *= e_factor
                        e_to_ee *= e_factor
                        e_to_ic *= e_factor
                        e_to_sec *= e_factor
                        coll_ex += e_to_ex
                        coll_ion += e_to_ion
                        coll_ex_HeI += e_to_ex_HeI
                        coll_ion_HeI += e_to_ion_HeI
                        coll_ex_HeII += e_to_ex_HeII
                        coll_ion_HeII += e_to_ion_HeII
                        coll_loss += e_to_ee
                        IC_total += e_to_ic
                        photons_particles += photons_particles_total_add*e_factor
                        e_factor_int = np.floor(e_to_sec / secondary_energy)
                        print electrons[i]*6.24e11, e_factor_int, e_to_sec*6.24e11, secondary_energy*6.24e11
                        if e_factor_int < 1:
                            e_factor_int = 1.0
                        # print electrons[i]*6.24e11, energy_to_distribute*6.24e11
                        electrons[i] -= energy_to_distribute
                        # print energy_to_distribute/(e_to_ex+e_to_ion+e_to_ex_HeI+e_to_ion_HeI+e_to_ex_HeII+e_to_ion_HeII+e_to_ic+e_to_ee)
                        # print energy_to_distribute/(e_to_ee+e_to_ex+e_to_ion+e_to_ic+e_to_sec+
                        #                                    e_to_ion_HeI+e_to_ion_HeII+
                        #                                    e_to_ex_HeI+e_to_ex_HeII)
                        secondary_electrons_list = rhoE_list[np.random.randint(1000, size=e_factor_int)]
                        secondary_electrons_list *= e_to_sec/np.sum(secondary_electrons_list)
                        electrons = np.hstack([electrons, secondary_electrons_list])
                    else:
                        tau_IC = electrons[i] / total_probability
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
                        if electrons[i] > cs.const_HI_ex / 6.24e11:
                            if rand_a[0] < tau_ex*tau:
                                electrons[i] -= cs.const_HI_ex / 6.24e11
                                coll_ex += cs.const_HI_ex / 6.24e11
                        if electrons[i] > cs.const_HeI_ex / 6.24e11:
                            if rand_a[1] < tau_ex_HeI*tau:
                                electrons[i] -= cs.const_HeI_ex / 6.24e11
                                coll_ex_HeI += cs.const_HeI_ex / 6.24e11
                        if electrons[i] > cs.const_HeII_ex / 6.24e11:
                            if rand_a[2] < tau_ex_HeII*tau:
                                electrons[i] -= cs.const_HeII_ex / 6.24e11
                                coll_ex_HeII += cs.const_HeII_ex / 6.24e11
                        if electrons[i] > cs.const_HI_ion / 6.24e11:
                            if rand_a[3] < tau_ion*tau:
                                temp = cs.rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= cs.const_HI_ion / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion += cs.const_HI_ion / 6.24e11
                        if electrons[i] > cs.const_HeI_ion / 6.24e11:
                            if rand_a[4] < tau_ion_HeI*tau:
                                temp = cs.rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= cs.const_HeI_ion / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion_HeI += cs.const_HeI_ion / 6.24e11
                        if electrons[i] > cs.const_HeII_ion / 6.24e11:
                            # print tau_ion_HeII*tau, tau_ion*tau
                            if rand_a[5] < tau_ion_HeII*tau:
                                temp = cs.rhoE(np.array([electrons[i]]), 8.0/6.24e11)[0]
                                electrons[i] -= cs.const_HeII_ion / 6.24e11
                                electrons[i] = electrons[i]-temp
                                electrons = np.append(electrons, temp)
                                coll_ion_HeII += cs.const_HeII_ion / 6.24e11
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
            inter_mask = np.where((electrons < E0_list[i_E-1]/6.24e11) & (electrons < E_extr_thresh/6.24e11))[0]# & (electrons < Eth_IC/6.24e11))[0]
            IC_total_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['IC'])
            # IC_soft_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 4])
            # IC_hard_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results[:, 5])
            coll_ion_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['H_I_ion'])
            coll_ex_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['H_I_ex'])
            coll_ion_HeI_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['He_I_ion'])
            coll_ex_HeI_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['He_I_ex'])
            coll_ion_HeII_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['He_II_ion'])
            coll_ex_HeII_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['He_II_ex'])
            coll_loss_frac = np.interp(electrons[inter_mask], E0_list/6.24e11, results['ee'])
            IC_total += np.sum(IC_total_frac * electrons[inter_mask])
            # IC_soft += np.sum(IC_soft_frac * electrons[inter_mask])
            # IC_hard += np.sum(IC_hard_frac * electrons[inter_mask])
            coll_ion += np.sum(coll_ion_frac * electrons[inter_mask])
            coll_ex += np.sum(coll_ex_frac * electrons[inter_mask])
            coll_ion_HeI += np.sum(coll_ion_HeI_frac * electrons[inter_mask])
            coll_ex_HeI += np.sum(coll_ex_HeI_frac * electrons[inter_mask])
            coll_ion_HeII += np.sum(coll_ion_HeII_frac * electrons[inter_mask])
            coll_ex_HeII += np.sum(coll_ex_HeII_frac * electrons[inter_mask])
            coll_loss += np.sum(coll_loss_frac * electrons[inter_mask])
            electrons = np.delete(electrons, inter_mask)
    results['IC'][i_E] = IC_total*6.24e11 / (MC_N*E0_list[i_E])
    results['H_I_ion'][i_E] = coll_ion*6.24e11 / (MC_N*E0_list[i_E])
    results['H_I_ex'][i_E] = coll_ex*6.24e11 / (MC_N*E0_list[i_E])
    results['He_I_ion'][i_E] = coll_ion_HeI*6.24e11 / (MC_N*E0_list[i_E])
    results['He_I_ex'][i_E] = coll_ex_HeI*6.24e11 / (MC_N*E0_list[i_E])
    results['He_II_ion'][i_E] = coll_ion_HeII*6.24e11 / (MC_N*E0_list[i_E])
    results['He_II_ex'][i_E] = coll_ex_HeII*6.24e11 / (MC_N*E0_list[i_E])
    results['ee'][i_E] = coll_loss*6.24e11 / (MC_N*E0_list[i_E])
    if results['IC'][i_E]/E0 > 0.9:
        E_extr_thresh = E0
    photons_particles_all[:, i_E] = photons_particles / MC_N
    print np.sum(photons_particles_all[:, i_E]*EgeV_list)
    print E0, E_extr_thresh #, np.sum((results[i_E]))


np.savez('output/%05i-%01.5f-%08.1f-snap'%(z,xiH,delta) + '_'+cs.cs['collion']+'_'+cs.cs['collex']+'_'+cs.cs['photion'] +'.npz', E0_list=E0_list, results=results, photons_particles_all=photons_particles_all, Eg_list=Eg_list, EgeV_list=EgeV_list)

# plt.plot(E0_list, np.sum(results[:, 1:4], 1))
# plt.plot(E0_list, results['ee'],'--k')
plt.plot(E0_list, results['H_I_ion'],'b')
plt.plot(E0_list, results['H_I_ex'],'--b')
plt.plot(E0_list, results['He_I_ion'],'r')
plt.plot(E0_list, results['He_I_ex'],'--r')
plt.plot(E0_list, results['He_II_ion'], 'k')
plt.plot(E0_list, results['He_II_ex'], '--k')
plt.plot(E0_list, results['IC'], '-')
# plt.plot(E0_list, results[:, 1])
# plt.plot(E0_list, results[:, 2])
# plt.plot(E0_list, results[:, 3])
plt.xscale('log')
plt.yscale('log')
plt.plot(E0_list, results['H_I_ion']+results['H_I_ex']+
                  results['He_I_ion']+results['He_I_ex']+
                  results['He_II_ion']+results['He_II_ex']+results['ee'], '--k')
plt.plot(E0_list, results['IC'], '-')
plt.plot(E0_list, results['IC']+results['H_I_ion']+results['H_I_ex']+results['He_I_ion']+results['He_I_ex']+results['He_II_ion']+results['He_II_ex']+results['ee'], '--k')
plt.plot([1e1,1e2,1e3,1e4,1e5,1e6], [0,0,0,0.35,0.85,1.0],'o')
plt.xscale('log')
plt.xlabel(r'$E\;\mathrm{[eV]}$')
plt.show()

