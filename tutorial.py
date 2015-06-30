__author__ = 'A'

import numpy as np
import matplotlib.pyplot as plt

x = np.ones(100)*9.
y = np.arange(1, 100, 2)

r = np.linspace(5e5, 10e5, 100)

plt.plot(r, r**2, 'o')
plt.yscale('log')
plt.xlabel('abc')

C = np.ones(100)
K = np.ones(100)
T = np.ones(100)*10**8

def Capacity(T, rho):
    return T**2-rho

dt = 0.001

for i in np.arange(1, 100):
    C = Capacity(T, rho)

    print i





