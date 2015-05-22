__author__ = 'kaurov'

# Furl 2010
def t_cool(E,z,xi=0):
    return 5e5*(1-xi)**-1*(E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Furl 2010
def t_loss_e(E, z, xi=0):
    return 5e3*(xi+1e-10)**-1 * (E/1000.)**1.5*((1.+z)/10.)**-3*3.14e7

# Shull 1985
def sigmaHe(E):
    return 2.75e-15 * np.log(E/13.6*6.24e11) / (E*6.24e11)
# def omegap(n, m):
    # return np.sqrt(4.0*np.pi*n*)


Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11

electrons = np.zeros(1)
E0 = 1e3
electrons[0] = E0 / 6.24e11
# electrons /= np.trapz(electrons, Eg_list)
# print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)
z=10
T_CMB=2.7*(1+z)

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

tau = 3.14e7*1e5
T_tau=0

photons_particles = np.zeros(len(Eg_list))

coll_ion_N = 0
coll_ion = 0
coll_loss = 0
coll_ex = 0
# Do electrons-photons interactions
for j in range(100):
    chances = np.random.random([len(electrons), 2])
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
        xi = 0.01
        tau_cool = t_cool(electrons[i]*6.24e11,z,xi)
        tau_loss_e = t_loss_e(electrons[i]*6.24e11,z,xi)
        tau_ex = 1./((1-xi)*nb)/sigmaHe(electrons[i])/c
        if total_probability > 0:
            tau = electrons[i]/total_probability/25.0
        else:
            total_probability = 0
            tau = 0
            if electrons[i] > 13.6/6.24e11:
                tau = np.min(np.array([tau_cool, tau_loss_e, tau_ex]))/25.0
            else:
                tau = np.min(np.array([tau_loss_e, tau_ex]))/25.0
        T_tau += tau
        electrons[i] -= total_probability*tau
        coll_ion_N += 1.0*tau/tau_cool
        if electrons[i] > 13.6/6.24e11:
            coll_ion += 1.0*tau/tau_cool*electrons[i]
            electrons[i] -= tau/tau_cool*electrons[i]
        coll_loss += 1.0*tau/tau_loss_e*electrons[i]
        electrons[i] -= tau/tau_loss_e*electrons[i]
        coll_ex += 1.0*tau/tau_ex*electrons[i]
        electrons[i] -= tau/tau_ex*electrons[i]
        photons_particles += photons_particles_total_add*tau
        print coll_ion / (coll_ion + coll_loss + coll_ex)
        print coll_loss / (coll_ion + coll_loss + coll_ex)
        print coll_ex / (coll_ion + coll_loss + coll_ex)
    # if j % 200 == 199:
        # plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)


print np.log10(np.sum(photons_particles*Eg_list*6.24e11))
print np.log10(coll_ion*6.24e11)
print np.log10(coll_loss*6.24e11)
print np.log10(np.sum(photons_particles*Eg_list*6.24e11)+electrons[0]*6.24e11+(coll_ion+coll_loss)*6.24e11)

print coll_ion / (coll_ion + coll_loss + coll_ex)

