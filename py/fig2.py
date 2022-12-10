#!/usr/bin/env python3
#==============================================================================
# fig2.py
#
# Plot the output saved by 'results.py'.
#
# Author: Stanley A. Baronett
# Created: 2022-12-09
# Updated: 2022-12-09
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np

# Load and plot data
npz = np.load('../npz/results.npz')
ns, xs, zs, rhos = npz['ns'], npz['xs'], npz['zs'], npz['rhos']
fig, axs = plt.subplots(3, sharex=True, figsize=(4.33, 5.42))

axs[0].plot(ns, xs)
axs[1].semilogy(ns, -zs)
axs[2].semilogy(ns, rhos)

# Format and save figure
for ax in axs.flat:
    ax.grid()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True)

axs[0].set(ylabel=r'$\xi_1$')
axs[1].set(ylabel=r'$-\theta_n^\prime(\xi_1)$', ylim=(6e-3, 1.2))
axs[2].set(xlabel=r'$n$', ylabel=r'$\rho_\mathrm{c} / \langle\rho\rangle$',
           ylim=(8e-1, 1.2e3))
plt.subplots_adjust(hspace=0)
plt.savefig(f'../figs/fig2.pdf', bbox_inches='tight', pad_inches=0.01)
plt.show()
