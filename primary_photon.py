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

cs = cross_sections(cs={'photion': 'VFKY1996', 'collion': 'AR', 'collex': 'SKD'})

Eg_list = np.logspace(-2, 12, 1000)/6.24e11
temp = np.zeros([len(Eg_list), 4])
EgeV_list = Eg_list*6.24e11

from scipy.interpolate import InterpolatedUnivariateSpline

crossections = np.loadtxt('radiator/datafiles/crossdata.dat')
crossectionsHe = np.loadtxt('radiator/datafiles/crossdataHe.dat')
crossectionsO = np.loadtxt('radiator/datafiles/crossdataO.dat')
crossectionsHe2 = np.loadtxt('radiator/datafiles/photoionHe.dat')
s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,-1]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '-', lw=1, label='H cross total')
s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
factor = 10**s(np.log10(EgeV_list/1e6))
plt.plot(EgeV_list, factor*1e-24, '-', lw=1, label='O cross total')
plt.plot(EgeV_list, cs.sigmaX(EgeV_list, 1, 1)*1e-18, label='H phot 1996')
plt.plot(EgeV_list, cs.sigmaX(EgeV_list, 8, 8)*1e-18, label='O phot 1996')
plt.xscale('log')
plt.yscale('log')
plt.legend()

plt.plot(EgeV_list, cs.sigmaX(EgeV_list, 1, 1) / cs.sigmaX(EgeV_list, 8, 8))
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


Tau_list = np.zeros(100000)

nb = 2.2e-7
# Eg_list = np.logspace(-6, 12, 1000)/6.24e11
nu_list = Eg_list / (hbar*2*np.pi)
EgeV_list = Eg_list*6.24e11

# Sigma_pair_prod = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,4]+crossections[:,5]))*1e-24
# Sigma_pair_prod[np.isnan(Sigma_pair_prod)] = 0.0
#
# # s = InterpolatedUnivariateSpline(np.log10(crossections[:,0]), np.log10(crossections[:,3]), k=1)
# # Sigma_photoion = 10**s(np.log10(EgeV_list/1e6))*1e-24
# Sigma_photoion = sigmaX(EgeV_list, 1, 1)*1e-18
# Sigma_photoion[np.isnan(Sigma_photoion)] = 0.0
# Sigma_photoion_He = sigmaX(EgeV_list, 2, 2)*1e-18
# Sigma_photoion_He[np.isnan(Sigma_photoion_He)] = 0.0
# Sigma_photoion_O = sigmaX(EgeV_list, 8, 8)*1e-18
#
# plt.plot(EgeV_list, Sigma_photoion)
# plt.plot(EgeV_list, Sigma_photoion_He)
# plt.plot(EgeV_list, Sigma_photoion_O)
#
# plt.xscale('log'); plt.yscale('log')
#
# Sigma_collisional_ion = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,2]))*1e-24
# Sigma_collisional_ion[np.isnan(Sigma_collisional_ion)] = 0.0
#
# # s = InterpolatedUnivariateSpline(np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,-1]), k=1)
# # Sigma_photoion_He = 10**s(np.log10(EgeV_list/1e6))*1e-24
# # Sigma_photoion_He[np.isnan(Sigma_photoion_He)] = 0.0
#
# Sigma_collisional_ion_He = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,2]))*1e-24
# Sigma_collisional_ion_He[np.isnan(Sigma_collisional_ion_He)] = 0.0
#
# Sigma_collisional_scatter = 10**np.interp(np.log10(EgeV_list/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,1]))*1e-24
# Sigma_collisional_scatter[np.isnan(Sigma_collisional_scatter)] = 0.0
#
# # s = InterpolatedUnivariateSpline(np.log10(crossectionsO[:,0]), np.log10(crossectionsO[:,-1]), k=1)
# # Sigma_photoion_O = 10**s(np.log10(EgeV_list/1e6))*1e-24
# # Sigma_photoion_O[np.isnan(Sigma_photoion_He)] = 0.0



z_end = 0.0
import cosmolopy as cp

tau = 3.14e7*1e0
photons_particles = np.zeros(len(Eg_list))
temp = 0


def runsim(E_ph_initial, z_start, bins=100):
    z_end = 0.0
    z_list = np.logspace(np.log10(z_start+1), np.log10(z_end+1), bins)-1.0
    age_list = cp.distance.age(z_list, **cp.fidcosmo)
    results = np.zeros(len(z_list))
    resultsH = np.zeros(len(z_list))
    resultsHe = np.zeros(len(z_list))
    total_mocks = 1000
    E_redshift_total = np.zeros(bins)
    E_other_total = np.zeros(bins)
    E_other_total_H = np.zeros(bins)
    E_other_total_He = np.zeros(bins)
    E_remaining_total = np.zeros(bins)
    for iii in range(total_mocks):
        randomii = np.random.random([len(z_list), 10])
        E_ph = 1.0*E_ph_initial
        E_redshift = 0.
        E_other = 0.
        ii=0
        while (E_ph>0) and (ii < (len(z_list)-2)):
            ii = ii+1
            E_ph_0 = 1.0*E_ph
            z = z_list[ii]
            tau = age_list[ii+1] - age_list[ii]
            dt = tau
            # column density
            N_naked = nb*c*tau*(1+z)**3
            N_H = N_naked*0.76
            N_He = N_naked*0.23/4
            N_O = N_naked*0.01/8
            # pair production
            temp_factor = (10**np.interp(np.log10(E_ph/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,4]+crossections[:,5]))*1e-24*N_naked)
            if np.isnan(temp_factor):
                chance_pp = 0
            else:
                chance_pp = 1.0 - np.exp(-temp_factor)

            # photoion
            opt_dep_photoion_H = (cs.sigmaX([E_ph], 1, 1)[0]*1e-18*N_H)
            opt_dep_photoion_He = (cs.sigmaX([E_ph], 2, 2)[0]*1e-18*N_He)
            opt_dep_photoion_O = (cs.sigmaX([E_ph], 8, 8)[0]*1e-18*N_O)
            chance_photoion_H = 1.0 - np.exp(-opt_dep_photoion_H)
            chance_photoion_He = 1.0 - np.exp(-opt_dep_photoion_He)
            chance_photoion_O = 1.0 - np.exp(-opt_dep_photoion_O)

            # collisional ion
            temp_factor_H = (10**np.interp(np.log10(E_ph/1e6), np.log10(crossections[:,0]), np.log10(crossections[:,2]))*1e-24*N_H)
            chance_coll_ion_H = 1.0 - np.exp(-temp_factor_H)
            temp_factor_He = (10**np.interp(np.log10(E_ph/1e6), np.log10(crossectionsHe[:,0]), np.log10(crossectionsHe[:,2]))*1e-24*N_H)
            chance_coll_ion_He = 1.0 - np.exp(-temp_factor_He)

            E_other = 0
            E_other_H = 0
            E_other_He = 0

            if randomii[ii, 0] < chance_pp+chance_photoion_H+chance_coll_ion_H+chance_photoion_He+chance_coll_ion_He:
                # print iii, ii
                tempH = chance_pp+chance_photoion_H+chance_coll_ion_H
                tempHe = chance_photoion_He+chance_coll_ion_He
                tempHHe = tempH + tempHe
                tempH /= tempHHe
                tempHe /= tempHHe
                resultsH[ii] += tempH
                resultsHe[ii] += tempHe
                results[ii] += 1
                E_other_H = tempH * E_ph
                E_other_He = tempHe * E_ph
                E_other = 1.0 * E_ph
                E_ph = 0

            # redshift
            E_ph = E_ph * (z_list[ii+1]+1) / (z_list[ii]+1)
            # print np.log10(E_ph),
            E_redshift = E_ph_0 - E_ph

            E_redshift_total[ii] += E_redshift
            E_remaining_total[ii] += E_ph
            E_other_total[ii] += E_other
            E_other_total_H[ii] += E_other_H
            E_other_total_He[ii] += E_other_He
    return results, E_redshift_total, E_remaining_total, E_other_total, E_other_total_H, E_other_total_He


bins=1000
total_mocks=1000

results430 = runsim(1e4, 50, bins=bins)
results630 = runsim(1e6, 50, bins=bins)
results930 = runsim(1e9, 50, bins=bins)

results41100 = runsim(1e4, 1100, bins=bins)
results61100 = runsim(1e6, 1100, bins=bins)
results91100 = runsim(1e9, 1100, bins=bins)

print E_redshift_total, E_remaining_total, E_other_total



f,(ax,ax2) = plt.subplots(1,2,sharey=True)

z_list = np.logspace(np.log10(50+1), np.log10(z_end+1), bins)-1.0
ax.plot(z_list, np.cumsum(results430[3])/total_mocks/1e4, 'k:', lw=1)
ax.plot(z_list, np.cumsum(results630[3])/total_mocks/1e6, 'k--', lw=1)
ax.plot(z_list, np.cumsum(results930[3])/total_mocks/1e9, 'k-', lw=1)
z_list = np.logspace(np.log10(1100+1), np.log10(z_end+1), bins)-1.0
ax.plot(z_list, np.cumsum(results41100[3])/total_mocks/1e4, 'k:', lw=1)
ax.plot(z_list, np.cumsum(results61100[3])/total_mocks/1e6, 'k--', lw=1)
ax.plot(z_list, np.cumsum(results91100[3])/total_mocks/1e9, 'k-', lw=1)


z_list = np.logspace(np.log10(50+1), np.log10(z_end+1), bins)-1.0
ax2.plot(z_list, np.cumsum(results430[3])/total_mocks/1e4, 'k:', lw=1)
ax2.plot(z_list, np.cumsum(results630[3])/total_mocks/1e6, 'k--', lw=1)
ax2.plot(z_list, np.cumsum(results930[3])/total_mocks/1e9, 'k-', lw=1)
z_list = np.logspace(np.log10(1100+1), np.log10(z_end+1), bins)-1.0
ax2.plot(z_list, np.cumsum(results41100[3])/total_mocks/1e4, 'k:', lw=1)
ax2.plot(z_list, np.cumsum(results61100[3])/total_mocks/1e6, 'k--', lw=1)
ax2.plot(z_list, np.cumsum(results91100[3])/total_mocks/1e9, 'k-', lw=1)



ax.set_xlim(0,54) # outliers only
ax2.set_xlim(100,1150) # most of the dat# a
ax.set_ylim(-0.01,1.01) # outliers only
ax2.set_ylim(-0.01,1.01) # outliers only

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
# ax.xaxis.tick_top()
ax2.tick_params(labeltop='off') # don't put tick labels at the top

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)      # top-left diagonal
ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
ax2.plot((-d,+d),(1-d,1+d), **kwargs)

ax.set_ylabel(r'$\mathrm{Absorbed\; energy\; fraction}$')
ax.set_xlabel(r'$\mathrm{Redshift}$')

plt.subplots_adjust(wspace=0.1)

plt.savefig('photon_ans.pdf')

#######################################

results430 = runsim(1e6, 30, bins=bins)

f,(ax,ax2) = plt.subplots(1,2,sharey=True)

z_list = np.logspace(np.log10(30+1), np.log10(z_end+1), bins)-1.0
ax.plot(z_list, np.cumsum(results430[3])/total_mocks/1e4, 'k:', lw=1)
ax.plot(z_list, np.cumsum(results430[4])/total_mocks/1e4, 'k--', lw=1)
ax.plot(z_list, np.cumsum(results430[5])/total_mocks/1e4, 'k-', lw=1)
# z_list = np.logspace(np.log10(1100+1), np.log10(z_end+1), bins)-1.0
# ax.plot(z_list, np.cumsum(results41100[3])/total_mocks/1e4, 'k:', lw=1)
# ax.plot(z_list, np.cumsum(results41100[4])/total_mocks/1e4, 'k--', lw=1)
# ax.plot(z_list, np.cumsum(results41100[5])/total_mocks/1e4, 'k-', lw=1)


z_list = np.logspace(np.log10(30+1), np.log10(z_end+1), bins)-1.0
ax2.plot(z_list, np.cumsum(results430[3])/total_mocks/1e4, 'k:', lw=1)
ax2.plot(z_list, np.cumsum(results430[4])/total_mocks/1e4, 'k--', lw=1)
ax2.plot(z_list, np.cumsum(results430[5])/total_mocks/1e4, 'k-', lw=1)
# z_list = np.logspace(np.log10(1100+1), np.log10(z_end+1), bins)-1.0
# ax2.plot(z_list, np.cumsum(results41100[3])/total_mocks/1e4, 'k:', lw=1)
# ax2.plot(z_list, np.cumsum(results41100[4])/total_mocks/1e4, 'k--', lw=1)
# ax2.plot(z_list, np.cumsum(results41100[5])/total_mocks/1e4, 'k-', lw=1)



ax.set_xlim(0,34) # outliers only
ax2.set_xlim(100,1150) # most of the dat# a
ax.set_ylim(-0.01,1.01) # outliers only
ax2.set_ylim(-0.01,1.01) # outliers only

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
# ax.xaxis.tick_top()
ax2.tick_params(labeltop='off') # don't put tick labels at the top

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)      # top-left diagonal
ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d),(-d,+d), **kwargs)   # bottom-left diagonal
ax2.plot((-d,+d),(1-d,1+d), **kwargs)

ax.set_ylabel(r'$\mathrm{Fraction}$')
ax.set_xlabel(r'$\mathrm{Redshift}$')

plt.subplots_adjust(wspace=0.1)

#######################################

bins=1000
total_mocks=10000
z_list = np.logspace(np.log10(30+1), np.log10(z_end+1), bins)-1.0

plt.plot(z_list, np.cumsum(E_other_total930)/total_mocks/1e9, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total930)/total_mocks/1e9, 'k--', lw=1)

plt.plot(z_list, np.cumsum(E_other_total630)/total_mocks/1e6, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total630)/total_mocks/1e6, 'k--', lw=1)

plt.plot(z_list, np.cumsum(E_other_total430)/total_mocks/1e4, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total430)/total_mocks/1e4, 'k--', lw=1)




bins=1000
total_mocks=10000
z_list = np.logspace(np.log10(1100+1), np.log10(z_end+1), bins)-1.0


plt.plot(z_list, np.cumsum(E_other_total91100)/total_mocks/1e9, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total91100)/total_mocks/1e9, 'k--', lw=1)

plt.plot(z_list, np.cumsum(E_other_total71100)/total_mocks/1e7, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total71100)/total_mocks/1e7, 'k--', lw=1)

plt.plot(z_list, np.cumsum(E_other_total61100)/total_mocks/1e6, 'k-', lw=1)
plt.plot(z_list, np.cumsum(E_redshift_total61100)/total_mocks/1e6, 'k--', lw=1)




###


bins=100
total_mocks=10000

E_ph_interaction_list = np.concatenate([np.logspace(np.log10(13.61), 2, 100)[:-1], np.logspace(2, 12, 30)])
Interaction_fraction = E_ph_interaction_list.copy() * 0.0
Interaction_fraction_H = E_ph_interaction_list.copy() * 0.0
Interaction_fraction_He = E_ph_interaction_list.copy() * 0.0
for i in range(len(E_ph_interaction_list)):
    print i
    results350 = runsim(E_ph_interaction_list[i], 30, bins=bins)
    Interaction_fraction[i] = np.sum(results350[3])/total_mocks/E_ph_interaction_list[i]
    Interaction_fraction_H[i] = np.sum(results350[4])/total_mocks/E_ph_interaction_list[i]
    Interaction_fraction_He[i] = np.sum(results350[5])/total_mocks/E_ph_interaction_list[i]

np.savez('E_effective.npz', E_ph_interaction_list = E_ph_interaction_list, Interaction_fraction=Interaction_fraction
         , Interaction_fraction_H=Interaction_fraction_H
         , Interaction_fraction_He=Interaction_fraction_He)

plt.plot(E_ph_interaction_list, Interaction_fraction_H/13.6)
plt.plot(E_ph_interaction_list, Interaction_fraction_He/24.4)
plt.xscale('log')
