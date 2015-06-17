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

def sigmaee(E, ne):
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

def eedEdt(E, ne):
    omega = c*np.sqrt(1-511e3**2/(511e3+E*6.24e11)**2)
    lnL = np.log(4*E*6.24e11/7.40e-11*ne)
    return 4*np.pi*ne*1.6e-20**4/9.1e-28/omega*lnL#*1e35
# def omegap(n, m):
    # return np.sqrt(4.0*np.pi*n*)


z=10
T_CMB=2.7*(1+z)

xi = 0.01


Eg_list = np.logspace(-6, 12, 100)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11


CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)


electrons = np.zeros(1)

E0_list = 13.6*10**np.linspace(0, 8, 100)

tau_list = np.zeros([len(E0_list), 10000]) * np.nan
IC_total_list = np.zeros([len(E0_list), 10000])
coll_ion_list = np.zeros([len(E0_list), 10000])
coll_ex_list = np.zeros([len(E0_list), 10000])
coll_loss_list = np.zeros([len(E0_list), 10000])

for iii in range(len(E0_list)):
    E0 = E0_list[iii]
    electrons[0] = E0 / 6.24e11
    tau = 3.14e7*1e-5
    T_tau=0
    photons_particles = np.zeros(len(Eg_list))
    coll_ion_N = 0
    coll_ion = 0
    coll_loss = 0
    coll_ex = 0
    IC_total = 0
    # Do electrons-photons interactions
    j = 0
    while (j < 9999) and (electrons[0] > 0.01 * E0 / 6.24e11) and (electrons[0] > 13.7/6.24e11):
        j = j+1
        # chances = np.random.random([len(electrons), 2])
        for i in range(len(electrons)):
            coll_ion = 0
            coll_loss = 0
            coll_ex = 0
            IC_total = 0
            gamma = electrons[i] / (me*c**2)
            sigma_IC = np.zeros(len(Eg_list))
            total_probability = 0
            photons_particles_total_add = 0
            for j1 in range(len(Eg_list)):
                sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                # sigma_IC[j1] = np.trapz(sigma_temp, Eg_list)
                # sigma_IC_total = sigma_IC[j1]*CMBphotons[j1]*np.gradient(Eg_list)[j1]
                # probability = sigma_IC_total * c * tau
                temp = sigma_temp * c * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                # print np.sum(temp),
                # print probability
                # if probability<=0:
                #     probability = 1e100
                # tau = electrons[i]/probability/10.0
                # tau=1e8
                probability = np.sum(temp*Eg_list)
                total_probability += probability
                photons_particles_total_add += temp
                # print probability
                # sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
                # sigma_temp = np.cumsum(sigma_temp)
                # sigma_temp /= sigma_temp.max()
                # # print j2
                # j2 = np.where(chances[i, 1] > sigma_temp)[0][-1]
                # photons_particles[j2] += 1
                # electrons[i] -= Eg_list[j2]
            tau_coll = t_coll(electrons[i]*6.24e11,z,xi)
            eedEdt_now = eedEdt(electrons[i], xi*nb*(1+z)**3)
            # tau_loss_e = t_loss_e(electrons[i]*6.24e11,z,xi)
            # tau_loss_e = electrons[i] / eedEdt_now
            tau_ex = 1./((1-xi)*nb)/sigmaHe(electrons[i])/(c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2))
            tau_loss_e = 1./((xi)*nb)/sigmaee(electrons[i], xi*nb*(1+z)**3)/(c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2))
            if total_probability > 0:
                tau_IC = electrons[i]/total_probability
            else:
                total_probability = 0
                tau_IC = 3.14e100
            if electrons[i] > 13.6/6.24e11:
                tau = np.min(np.array([tau_loss_e/0.05, tau_ex*electrons[i]*6.24e11/13.6, tau_IC]))/20.0
            else:
                tau = np.min(np.array([tau_loss_e/0.05, tau_ex*electrons[i]*6.24e11/13.6, tau_IC]))/20.0
            T_tau += tau
            e0z = electrons[i].copy()
            electrons[i] -= total_probability*tau
            IC_total += total_probability*tau
            if electrons[i] > 13.6/6.24e11:
                coll_ion += 1.0*tau/tau_ex*13.6/6.24e11
                electrons[i] -= tau/tau_ex*13.6/6.24e11
            if electrons[i] > 10.6/6.24e11:
                coll_ex += 1.0*tau/tau_ex*10.6/6.24e11
                electrons[i] -= tau/tau_ex*10.6/6.24e11
            #coll_ex += 1.0*tau/tau_ex*10.6/6.24e11
            #electrons[i] -= tau/tau_ex*10.6/6.24e11
            coll_loss += 1.0*tau/tau_loss_e*electrons[i]*0.05
            electrons[i] -= tau/tau_loss_e*electrons[i]*0.05
            # coll_loss += 1.0*tau*eedEdt_now
            # electrons[i] -= tau*eedEdt_now
            # coll_ex += 1.0*tau/tau_ex*e0z
            # electrons[i] -= tau/tau_ex*e0z
            photons_particles += photons_particles_total_add*tau
            tau_list[iii, j] = T_tau/3.14e7
            print iii, j, electrons[i]/E0*6.24e11
            print IC_total, coll_ion, coll_ex, coll_loss
            IC_total_list[iii, j] = IC_total
            coll_ion_list[iii, j] = coll_ion
            coll_ex_list[iii, j] = coll_ex
            coll_loss_list[iii, j] = coll_loss
    j += 1
    IC_total_list[iii, j] = 0
    coll_ion_list[iii, j] = 0
    coll_ex_list[iii, j] = 0
    coll_loss_list[iii, j] = electrons[i]

        # if j % 200 == 199:
            # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)


total_energy = IC_total_list+coll_ion_list+coll_ex_list+coll_loss_list
plt.plot(E0_list, IC_total_list.sum(1)/total_energy.sum(1), '-k')
plt.plot(E0_list, coll_ion_list.sum(1)/total_energy.sum(1), '--k')
plt.plot(E0_list, coll_ex_list.sum(1)/total_energy.sum(1), '-.k')
plt.plot(E0_list, coll_loss_list.sum(1)/total_energy.sum(1), ':k')
plt.xscale('log')
plt.xlim([10,1e9])
# plt.yscale('log')
plt.xlabel(r'$E_0\;\mathrm{[eV]}$')
plt.ylabel(r'$\mathrm{Fractional\;contribution}$')
fixlogax(plt.gca(), 'x')
plt.savefig('figures/FracVsE0z10xi001.pdf')

############################

for iii in [99]:
    E0 = E0_list[iii]
    electrons[0] = E0 / 6.24e11
    tau = 3.14e7*1e-5
    T_tau=0
    photons_particles = np.zeros(len(Eg_list))
    coll_ion_N = 0
    coll_ion = 0
    coll_loss = 0
    coll_ex = 0
    IC_total = 0
    # Do electrons-photons interactions
    j = 0
    while (j < 9999) and (electrons[0] > 13.7/6.24e11):
        j = j+1
        # chances = np.random.random([len(electrons), 2])
        for i in range(len(electrons)):
            coll_ion = 0
            coll_loss = 0
            coll_ex = 0
            IC_total = 0
            gamma = electrons[i] / (me*c**2)
            sigma_IC = np.zeros(len(Eg_list))
            total_probability = 0
            photons_particles_total_add = 0
            for j1 in range(len(Eg_list)):
                sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                # sigma_IC[j1] = np.trapz(sigma_temp, Eg_list)
                # sigma_IC_total = sigma_IC[j1]*CMBphotons[j1]*np.gradient(Eg_list)[j1]
                # probability = sigma_IC_total * c * tau
                temp = sigma_temp * c * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                # print np.sum(temp),
                # print probability
                # if probability<=0:
                #     probability = 1e100
                # tau = electrons[i]/probability/10.0
                # tau=1e8
                probability = np.sum(temp*Eg_list)
                total_probability += probability
                photons_particles_total_add += temp
                # print probability
                # sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
                # sigma_temp = np.cumsum(sigma_temp)
                # sigma_temp /= sigma_temp.max()
                # # print j2
                # j2 = np.where(chances[i, 1] > sigma_temp)[0][-1]
                # photons_particles[j2] += 1
                # electrons[i] -= Eg_list[j2]
            tau_coll = t_coll(electrons[i]*6.24e11,z,xi)
            eedEdt_now = eedEdt(electrons[i], xi*nb*(1+z)**3)
            # tau_loss_e = t_loss_e(electrons[i]*6.24e11,z,xi)
            # tau_loss_e = electrons[i] / eedEdt_now
            tau_ex = 1./((1-xi)*nb)/sigmaHe(electrons[i])/(c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2))
            tau_loss_e = 1./((xi)*nb)/sigmaee(electrons[i], xi*nb*(1+z)**3)/(c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2))
            if total_probability > 0:
                tau_IC = electrons[i]/total_probability
            else:
                total_probability = 0
                tau_IC = 3.14e100
            tau = np.min(np.array([tau_loss_e/0.05, tau_ex*electrons[i]*6.24e11/13.6, tau_IC]))
            temp = (tau_IC / (tau_loss_e/0.05))
            if (temp>1e-7) & (temp < 1e4):
                tau = tau/1500.
            else:
                tau = tau/300.
            T_tau += tau
            e0z = electrons[i].copy()
            electrons[i] -= total_probability*tau
            IC_total += total_probability*tau
            if electrons[i] > 13.6/6.24e11:
                coll_ion += 1.0*tau/tau_ex*13.6/6.24e11
                electrons[i] -= tau/tau_ex*13.6/6.24e11
            if electrons[i] > 10.6/6.24e11:
                coll_ex += 1.0*tau/tau_ex*10.6/6.24e11
                electrons[i] -= tau/tau_ex*10.6/6.24e11
            #coll_ex += 1.0*tau/tau_ex*10.6/6.24e11
            #electrons[i] -= tau/tau_ex*10.6/6.24e11
            coll_loss += 1.0*tau/tau_loss_e*electrons[i]*0.05
            electrons[i] -= tau/tau_loss_e*electrons[i]*0.05
            # coll_loss += 1.0*tau*eedEdt_now
            # electrons[i] -= tau*eedEdt_now
            # coll_ex += 1.0*tau/tau_ex*e0z
            # electrons[i] -= tau/tau_ex*e0z
            photons_particles += photons_particles_total_add*tau
            tau_list[iii, j] = T_tau/3.14e7
            print iii, j, electrons[i]/E0*6.24e11
            print IC_total, coll_ion, coll_ex, coll_loss
            IC_total_list[iii, j] = IC_total
            coll_ion_list[iii, j] = coll_ion
            coll_ex_list[iii, j] = coll_ex
            coll_loss_list[iii, j] = coll_loss
    j += 1
    IC_total_list[iii, j] = 0
    coll_ion_list[iii, j] = 0
    coll_ex_list[iii, j] = 0
    coll_loss_list[iii, j] = electrons[i]

        # if j % 200 == 199:
            # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)

nnn = 99
plt.clf()
total_energy = IC_total_list+coll_ion_list+coll_ex_list+coll_loss_list
ax1=plt.subplot(211)
plt.plot(tau_list[nnn, :]/3.14e7, (IC_total_list/total_energy)[nnn,:], '-k')
plt.plot(tau_list[nnn, :]/3.14e7, (coll_ion_list/total_energy)[nnn,:], '--k')
plt.plot(tau_list[nnn, :]/3.14e7, (coll_ex_list/total_energy)[nnn,:], '-.k')
plt.plot(tau_list[nnn, :]/3.14e7, (coll_loss_list/total_energy)[nnn,:], ':k')
# plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-2,3e4])
# plt.yscale('log')
plt.ylabel(r'$\mathrm{Fraction}$')
ax1.set_xticklabels([])
plt.subplot(212)
plt.plot(tau_list[61, :]/3.14e7, 1.0-np.cumsum(total_energy[nnn, :])/np.sum(total_energy[nnn, :]),'-k'); #plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-2,3e4])
plt.ylabel(r'$E(t)/E_0$')
plt.xlabel(r'$\mathrm{Time}\;\mathrm{[yrs]}$')
fixlogax(plt.gca(), 'x')
plt.savefig('figures/FracVsTimez10xi001.pdf')