import numpy as np
import matplotlib.pyplot as plt
from radiator.figurestyle import *
# from rstat.crossections import sigmaH
define_figure_style()


files = ['output\\00100-0.01000-001000.0-snap.npz',\
         #'output\\00010-0.01000-000000.0-snap.npz',\
         'output\\00100-0.01000-000000.0-snap.npz']

styles = ['-', '--', ':']
linews = [2, 1, 1]

for i in range(len(files)):
    data = np.load(files[i])
    E0_list = data['E0_list']
    results = data['results']
    photons_particles_all = data['photons_particles_all']
    plt.plot(E0_list, np.sum(results[:, 1:4],1), styles[i], color='k', lw=1)
    # plt.plot(E0_list, results[:, 0])
    plt.plot(E0_list, results[:, 4], styles[i], color='k', lw=1)
    plt.plot(E0_list, results[:, 5], styles[i], color='k', lw=1)
    plt.xscale('log')
    plt.xlabel(r'$E\;\mathrm{[eV]}$')