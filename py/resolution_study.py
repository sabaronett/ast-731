#!/usr/bin/env python3
#==============================================================================
# resolution_study.py
#
# Compute and save the absolute differences of the polytropic parameters ξ₁ and
# θₙ'(ξ₁) of the Lane--Emden equataion with their analytic solutions given by 
# Hansen, Kawaler & Trimble (2004, Table 7.1).
#
# Author: Stanley A. Baronett
# Created: 2022-12-08
# Updated: 2022-12-09
#==============================================================================
import numpy as np
import shooting_method as sm
import sys

ns = np.asarray([0, 1])                               # polytropic indices
h_min, h_max, num = 1e-6, 1e-2, 200
hs = np.geomspace(h_min, h_max, num=num)              # range of step sizes
analytics = [[np.sqrt(6), -np.sqrt(6)/3],             # analytic solutions
             [np.pi,      -1/np.pi     ]]
dxs = np.zeros((len(ns), num))                        # to store δξ₁
dzs = np.zeros((len(ns), num))                        # to store δθₙ'(ξ₁)

print("Computing and storing δξ₁ and δθₙ'(ξ₁) for...", flush=True)
for i, n in enumerate(ns):
    xs, zs = np.zeros(num), np.zeros(num)             # to store ξ₁ and θₙ'(ξ₁)

    print(f'  n = {n}...', flush=True)

    for j, h in enumerate(hs):
        x = 1e-16                                     # ξ
        y = sm.yfunc(x, n)                            # θₙ
        z = sm.zfunc(x, n)                            # θₙ' = dθₙ/dξ = dy/dx

        while y.real > 0:
            x, y, z = sm.rk4(n, x, y, z, h)
        
        xs[j], zs[j] = x, z.real
        sys.stdout.write(f'\r    {(j+1)/hs.size:.1%}')

    dxs[i][:] = np.abs(np.asarray(xs) - analytics[i][0])
    dzs[i][:] = np.abs(np.asarray(zs) - analytics[i][1])
    print('  Done.', flush=True)

print("Saving results for plotting in 'fig1.py'...", flush=True)
np.savez_compressed('../npz/resolution_study', ns=ns, hs=hs, dxs=dxs, dzs=dzs)
print('Finished.', flush=True)
