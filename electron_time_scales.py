__author__ = 'kaurov'
import numpy as np

# Furl 2010
def t_coll(E,z,xi=0):
    return 5e5*(1-xi)**-1*(E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Furl 2010
def t_loss_e(E, z, xi=0):
    return 5e3*(xi+1e-10)**-1 * (E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Shull 1985
def sigmaHe(E):
    return 2.75e-15 * np.log(E/13.6*6.24e11) / (E*6.24e11)

def eedEdt(E, ne):
    omega = c*np.sqrt(1-511e3**2/(511e3+E*6.24e11)**2)
    lnL = np.log(4*E*6.24e11/7.40e-11*ne)
    return 4.0*np.pi*ne*1.6e-19**4/9.1e-28/omega*lnL*1e35
# def omegap(n, m):
    # return np.sqrt(4.0*np.pi*n*)


Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11

electrons = np.zeros(1)
E0 = 30
electrons[0] = E0 / 6.24e11
# electrons /= np.trapz(electrons, Eg_list)
# print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)
z=10
T_CMB=2.7*(1+z)

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

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
while (j < 100) and (electrons[0] > 13.7/6.24e11):
    j = j+1
    # chances = np.random.random([len(electrons), 2])
    for i in range(len(electrons)):
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
        xi = 0.1
        tau_coll = t_coll(electrons[i]*6.24e11,z,xi)
        eedEdt_now = eedEdt(electrons[i], xi*nb*(1+z)**3)
        # tau_loss_e = t_loss_e(electrons[i]*6.24e11,z,xi)
        tau_loss_e = electrons[i] / eedEdt_now
        tau_ex = 1./((1-xi)*nb)/sigmaHe(electrons[i])/(c*np.sqrt(1-511e3**2/(511e3+electrons[i]*6.24e11)**2))
        if total_probability > 0:
            tau_IC = electrons[i]/total_probability
        else:
            total_probability = 0
            tau_IC = 3.14e100
        if electrons[i] > 13.6/6.24e11:
            tau = np.min(np.array([tau_loss_e, tau_ex, tau_IC]))/100.0
        else:
            tau = np.min(np.array([tau_loss_e, tau_ex, tau_IC]))/100.0
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
        coll_loss += 1.0*tau*eedEdt_now
        electrons[i] -= tau*eedEdt_now
        # coll_ex += 1.0*tau/tau_ex*e0z
        # electrons[i] -= tau/tau_ex*e0z
        photons_particles += photons_particles_total_add*tau
        print T_tau/3.14e7
        print electrons[i]/E0*6.24e11
        print IC_total / (coll_ion + coll_loss + coll_ex + IC_total),
        print coll_ion / (coll_ion + coll_loss + coll_ex + IC_total),
        print coll_ex / (coll_ion + coll_loss + coll_ex + IC_total),
        print coll_loss / (coll_ion + coll_loss + coll_ex + IC_total)
    # if j % 200 == 199:
        # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)


print np.log10(np.sum(photons_particles*Eg_list*6.24e11))
print np.log10(coll_ion*6.24e11)
print np.log10(coll_loss*6.24e11)
print np.log10(np.sum(photons_particles*Eg_list*6.24e11)+electrons[0]*6.24e11+(coll_ion+coll_loss)*6.24e11)

print coll_ion / (coll_ion + coll_loss + coll_ex)

