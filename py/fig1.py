#!/usr/bin/env python3
#==============================================================================
# fig1.py
#
# Plot the output saved by 'resolution_study.py'.
#
# Author: Stanley A. Baronett
# Created: 2022-12-09
# Updated: 2022-12-09
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np

# Load and plot data
npz = np.load('../npz/resolution_study.npz')
ns, hs, dxs, dzs = npz['ns'], npz['hs'], npz['dxs'], npz['dzs']
fig, axs = plt.subplots(2, ns.size, sharex=True, sharey='row',
                        figsize=(8.66, 5.06))

for i in range(ns.size):
    axs[0][i].loglog(hs, dxs[i])
    axs[1][i].loglog(hs, dzs[i])

# Format and save figure
for ax in axs.flat:
    ax.grid()
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', top=True)

axs[0][0].set(title=r'$n=0$', ylabel=r'$\delta\xi_1$')
axs[0][1].set(title=r'$n=1$')
axs[1][0].set(xlabel=r'$h$', ylabel=r'$\delta\theta^\prime_n(\xi_1)$')
axs[1][1].set(xlabel=r'$h$')
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig(f'../figs/fig1.pdf', bbox_inches='tight', pad_inches=0.01)
plt.show()
