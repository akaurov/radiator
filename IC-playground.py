import numpy as np
import matplotlib.pyplot as plt
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

sigmaknV = np.vectorize(sigmakn, excluded=['Eg', 'gamma'])


Eg_list = np.logspace(-6, 12, 300)/6.24e11
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


electrons = np.zeros(1)
E00 = 1e9
electrons[0] = E00 / 6.24e11
# electrons /= np.trapz(electrons, Eg_list)
# print 'number of electrons per cm^3: ', np.trapz(electrons, Eg_list)
z=10
T_CMB=2.7*(1+z)

CMBphotons = 8*np.pi*2*np.pi*hbar*nu_list**3/c**3/(np.exp(2.0*np.pi*hbar*nu_list/kB/T_CMB)-1.0)/Eg_list/(hbar*2*np.pi)
print 'number of photons per cm^3: ', np.trapz(CMBphotons, Eg_list)

tau = 3.14e7*1e1
photons_particles = np.zeros(len(Eg_list))
photons_particles_add = np.zeros(len(Eg_list))

# Do electrons-photons interactions
for j in range(1000):
    chances = np.random.random([len(electrons), 2])
    for i in range(len(electrons)):
        if electrons[i]>10**5.42 / 6.24e11:
            gamma = electrons[i] / (me*c**2)
            sigma_IC = np.zeros(len(Eg_list))
            photons_particles_add = np.zeros(len(Eg_list))
            electrons_minus = 0
            for j1 in range(len(photons)):
                sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma) * np.gradient(Eg_list)
                # sigma_IC[j1] = np.trapz(sigma_temp, Eg_list)
                # sigma_IC_total = sigma_IC[j1]*CMBphotons[j1]*np.gradient(Eg_list)[j1]
                # probability = sigma_IC_total * c * tau
                temp = sigma_temp * c * CMBphotons[j1] * np.gradient(Eg_list)[j1]
                photons_particles_add += temp
                # print np.sum(temp),
                probability = np.sum(temp*Eg_list)
                electrons_minus += probability
                # print probability
                # sigma_temp = sigmakn(Eg_list, Eg_list[j1], gamma)
                # sigma_temp = np.cumsum(sigma_temp)
                # sigma_temp /= sigma_temp.max()
                # # print j2
                # j2 = np.where(chances[i, 1] > sigma_temp)[0][-1]
                # photons_particles[j2] += 1
                # electrons[i] -= Eg_list[j2]
            tau = electrons[i]/electrons_minus/100.0
            photons_particles_add *= tau
            photons_particles += photons_particles_add
            electrons_minus *= tau
            electrons[i] -= electrons_minus
            print tau, np.log10(electrons[0]*6.24e11)
    # if j % 200 == 199:
    #     plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)


print np.log10(np.sum(photons_particles*Eg_list*6.24e11))
print np.log10(np.sum(photons_particles*Eg_list*6.24e11)+electrons[0]*6.24e11)
plt.plot(EgeV_list, photons_particles/np.gradient(Eg_list)*Eg_list**2)
plt.plot(EgeV_list, CMBphotons*Eg_list**2)
plt.xscale('log')
plt.yscale('log')


def ICon30CMB(EgeV_list, E0, z):
    return 4e-2 * (EgeV_list/1e9)**0.5 * \
           np.exp(\
               -(EgeV_list/1e9)**1.25*5.4e3* \
               (1e9/E0)**2 / \
               ((1.+z)/61)**0.5 \
               )

plt.plot(EgeV_list, ICon30CMB(EgeV_list, E00, z),'--')
plt.xscale('log')
plt.yscale('log')