__author__ = 'A'
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

nb = 2.2e-7
# Furl 2010
def t_coll(E,z,xi=0):
    return 5e5*(1-xi)**-1*(E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Furl 2010
def t_loss_e(E, z, xi=0):
    return 5e3*(xi+1e-10)**-1 * (E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Shull 1985
def sigmaHe(E):
    return 2.75e-15 * np.log(E/13.6*6.24e11) / (E*6.24e11)

def sigmaee(E, ne, dE):
    f = dE/E
    omega = c*np.sqrt(1-511e3**2/(511e3+E*6.24e11)**2)
    lnL = np.log(4*E*6.24e11/7.40e-11*ne)
    return 0.3*7.82e-11*lnL *(E*6.24e11)**-2

### Sigma test

# E0_list = 13.6*10**np.linspace(0, 5, 100)
# plt.plot(E0_list, sigmaHe(E0_list/6.24e11))
# plt.plot(E0_list, sigmaee(E0_list/6.24e11, nb))
# plt.xscale('log')
# plt.yscale('log')

###

# Not complete!!!
def eedEdt(E, ne):
    omega = c*np.sqrt(1-511e3**2/(511e3+E*6.24e11)**2)
    lnL = 30.#np.log(4*E*6.24e11/7.40e-11*ne)
    return 4*np.pi*ne*4.8e-10**4/9.1e-28/omega*lnL#*1e35


z=100
T_CMB=2.7*(1+z)

xi = 0.004


Eg_list = np.logspace(-6, 12, 100)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11


CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)



# tau_list = np.zeros([len(E0_list), 10000]) * np.nan
# IC_total_list = np.zeros([len(E0_list), 10000])
# coll_ion_list = np.zeros([len(E0_list), 10000])
# coll_ex_list = np.zeros([len(E0_list), 10000])
# coll_loss_list = np.zeros([len(E0_list), 10000])

tau = 3.14e7*1e-5
T_tau=0
photons_particles = np.zeros(len(Eg_list))
coll_ion_N = 0
coll_ion = 0
coll_loss = 0
coll_ex = 0
IC_total = 0
# Do electrons-photons interactions
E0_list = np.logspace(1, 7, 30)
# E0_list = np.array([1e6])
results = np.zeros([len(E0_list), 10])
MC_N_list = np.zeros([len(E0_list)])

for i_E in range(len(E0_list)):
    E0 = E0_list[i_E]

    if E0<1e2:
        MC_N=1000
    elif E0<1e3:
        MC_N=100
    elif E0<1e5:
        MC_N=3
    else:
        MC_N=1
    MC_N_list[i_E]=MC_N

    coll_ion = 0
    coll_loss = 0
    coll_ex = 0
    IC_total = 0

    for iii in range(MC_N):
        electrons = np.zeros(1)
        electrons[0] = E0 / 6.24e11
        j=0
        while (electrons.sum()*6.24e11 > 1e-10):
            j = j+1
            # chances = np.random.random([len(electrons), 2])
            for i in range(len(electrons)):
                if electrons[i] > 10.2/6.24e11:
                    gamma = electrons[i] / (me*c**2)
                    v = c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2)
                    sigma_IC = np.zeros(len(Eg_list))
                    total_probability = 0
                    photons_particles_total_add = 0
                    if electrons[i] > 1e4/6.24e11:
                        for j1 in range(len(Eg_list)):
                            sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                            temp = sigma_temp * c * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                            probability = np.sum(temp*Eg_list)
                            total_probability += probability
                            photons_particles_total_add += temp
                    tau_coll = t_coll(electrons[i]*6.24e11, z, xi)
                    eedEdt_now = eedEdt(electrons[i], xi*nb*(1+z)**3)
                    tau_ex = (1-xi)*nb*(1+z)**3*sigmaHe(electrons[i])*v
                    # tau_loss_e = (xi)*nb*(1+z)**3*sigmaee(electrons[i], xi*nb*(1+z)**3)*v
                    # if total_probability > 0:
                    #     tau_IC = electrons[i]/total_probability
                    # else:
                    #     total_probability = 0
                    #     tau_IC = 3.14e100
                    if electrons[i] > 13.6/6.24e11:
                        tau = 0.1/tau_ex
                    T_tau += tau
                    e0z = electrons[i].copy()
                    electrons[i] -= total_probability*tau
                    IC_total += total_probability*tau

                    rand_a = np.random.rand(2)
                    if rand_a[0] < (1.0 - np.exp(-tau_ex*tau)):
                        electrons[i] -= 10.2 / 6.24e11
                        coll_ex += 10.2 / 6.24e11
                    if rand_a[1] < (1.0 - np.exp(-tau_ex*tau)):
                        electrons[i] -= 13.6 / 6.24e11
                        electrons[i] = electrons[i]/2.0
                        electrons = np.append(electrons, electrons[i])
                        coll_ion += 13.6 / 6.24e11
                    electrons[i] -= eedEdt_now * tau
                    coll_loss += eedEdt_now * tau

                    # photons_particles += photons_particles_total_add*tau
                    # tau_list[iii, j] = T_tau/3.14e7
                    # print i, j, electrons[i]/E0*6.24e11
                    # print IC_total*6.24e11, coll_ion*6.24e11, coll_ex*6.24e11, coll_loss*6.24e11
                else:
                    coll_loss += electrons[i]
                    electrons[i] = 0
            # print j

    results[i_E,:4] = IC_total*6.24e11, coll_ion*6.24e11, coll_ex*6.24e11, coll_loss*6.24e11
    print E0, results[i_E,:4]

plt.plot(E0_list, results[:, 0]/E0_list/MC_N_list)
plt.plot(E0_list, results[:, 1]/E0_list/MC_N_list)
plt.plot(E0_list, results[:, 2]/E0_list/MC_N_list)
plt.plot(E0_list, results[:, 3]/E0_list/MC_N_list)
plt.xscale('log')
# plt.yscale('log')


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