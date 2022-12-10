#!/usr/bin/env python3
#==============================================================================
# results.py
#
# Compute and save the polytropic parameters ξ₁, θₙ'(ξ₁), and ρ_c/⟨ρ⟩ for a
# range of polytropic indices n.
#
# Author: Stanley A. Baronett
# Created: 2022-12-08
# Updated: 2022-12-09
#==============================================================================
import numpy as np
import shooting_method as sm
import sys

h = 1e-5                                                 # step size
num = 50                                                 # no. of n samples
xs = np.zeros(num + 1)                                   # to store ξ₁
zs = np.zeros(num + 1)                                   # to store θₙ'(ξ₁)
rhos = np.zeros(num + 1)                                 # to store ρ_c/⟨ρ⟩
n_min, n_max = 0, 4                                      # range of n values
n_step = (n_max - n_min)/num
ns = [n_min]

print('Creating evenly-spaced range of n values...', flush=True)
for i in range(num):
    ns.append(ns[-1] + n_step)

print('  Done.\nComputing and storing parameters for each n...', flush=True)
for i, n in enumerate(ns):
    x = 1e-16                                            # ξ
    y = sm.yfunc(x, n)                                   # θₙ
    z = sm.zfunc(x, n)                                   # θₙ' = dθₙ/dξ = dy/dx

    while y.real > 0:
        x, y, z = sm.rk4(n, x, y, z, h)
    
    xs[i], zs[i], rhos[i] = x, z.real, -x/z.real/3
    sys.stdout.write(f'\r  {(i+1)/len(ns):.0%}')

print("  Done.\nSaving results for plotting in 'fig2.py'...", flush=True)
np.savez_compressed('../npz/results', ns=ns, xs=xs, zs=zs, rhos=rhos)

print("  Done.\nGenerating Table 1...", flush=True)
ns = [0, 1.0, 1.5, 2.0, 3.0, 4.0]

print()
print('--------------------------------')
print('Index n\t  ξ₁\t-θₙ′(ξ₁) ρ_c/<ρ>')
print('--------------------------------')

for n in ns:
    x = 1e-16
    y = sm.yfunc(x, n)
    z = sm.zfunc(x, n)

    while y.real > 0:
        x, y, z = sm.rk4(n, x, y, z, h)

    rhoc = -x/z.real/3

    print(f'  {n:.1f}\t{x:.4f}\t{-z.real:.5f}\t{rhoc:.4f}')

print()
print('Finished.', flush=True)
