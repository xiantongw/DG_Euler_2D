#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

es_bump0 = [1.73544377e-03, 1.97351555e-04, 1.38560815e-05]
es_bump1 = [1.12450757e-03, 3.96116861e-05, 2.79903085e-06]
es_bump2 = [6.43098298e-04, 8.51978470e-06, 1.86945460e-06]

cl_bump0 = [2.61354830e-01, 1.48070697e+00, 1.53622953e+00] 
cl_bump1 = [7.75166539e-01, 1.53101901e+00, 1.53741594e+00] 
cl_bump2 = [1.12214182e+00, 1.53759677e+00, 1.53695592e+00] 

cd_bump0 = [4.58931459e-01, 2.42078262e-02, 9.09358745e-04] 
cd_bump1 = [2.75198349e-01, 3.67324330e-03, 5.15474805e-05] 
cd_bump2 = [1.48905993e-01, 5.45915781e-04, 2.65279670e-05] 

cl_exact = 1.537095
cd_exact = 2.94278e-6
es_exact = 0.0

plt.figure(figsize=(28,5))
matplotlib.rcParams.update({'font.size': 16})

plt.subplot(1, 3, 1)
p = 0; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cl_p0 = np.abs(np.array([cl_bump0[0], cl_bump1[0], cl_bump2[0]]) - cl_exact)
r0 = np.polyfit(np.log(sqrt_dof), np.log(err_cl_p0), 1)[0]
plt.loglog(sqrt_dof, err_cl_p0, marker='s', label="p=0, k={:.2f}".format(r0))

p = 1; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cl_p1 = np.abs(np.array([cl_bump0[1], cl_bump1[1], cl_bump2[1]]) - cl_exact)
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_cl_p1), 1)[0]
plt.loglog(sqrt_dof, err_cl_p1, marker='s', label="p=1, k={:.2f}".format(r1))

p = 2; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cl_p2 = np.abs(np.array([cl_bump0[2], cl_bump1[2], cl_bump2[2]]) - cl_exact)
r2 = np.polyfit(np.log(sqrt_dof), np.log(err_cl_p2), 1)[0]
plt.loglog(sqrt_dof, err_cl_p2, marker='s',  label="p=2, k={:.2f}".format(r2))

plt.xlabel(r"$\sqrt{dof}$")
plt.ylabel(r"$Error\ of\ c_l$")
plt.legend(loc=3)

plt.subplot(1, 3, 2)
p = 0; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cd_p0 = np.abs(np.array([cd_bump0[0], cd_bump1[0], cd_bump2[0]]) - cd_exact)
r0 = np.polyfit(np.log(sqrt_dof), np.log(err_cd_p0), 1)[0]
plt.loglog(sqrt_dof, err_cd_p0, marker='s', label="p=0, k={:.2f}".format(r0))

p = 1; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cd_p1 = np.abs(np.array([cd_bump0[1], cd_bump1[1], cd_bump2[1]]) - cd_exact)
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_cd_p1), 1)[0]
plt.loglog(sqrt_dof, err_cd_p1, marker='s', label="p=1, k={:.2f}".format(r1))

p = 2; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_cd_p2 = np.abs(np.array([cd_bump0[2], cd_bump1[2], cd_bump2[2]]) - cd_exact)
r2 = np.polyfit(np.log(sqrt_dof[0:2]), np.log(err_cd_p2[0:2]), 1)[0]
plt.loglog(sqrt_dof, err_cd_p2, marker='s', label="p=2, k={:.2f}".format(r2))

plt.xlabel(r"$\sqrt{dof}$")
plt.ylabel(r"$Error\ of\ c_d$")
plt.legend(loc=3)


plt.subplot(1, 3, 3)
p = 0; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_es_p0 = np.abs(np.array([es_bump0[0], es_bump1[0], es_bump2[0]]) - es_exact)
r0 = np.polyfit(np.log(sqrt_dof), np.log(err_es_p0), 1)[0]
plt.loglog(sqrt_dof, err_es_p0, marker='s', label="p=0, k={:.2f}".format(r0))

p = 1; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_es_p1 = np.abs(np.array([es_bump0[1], es_bump1[1], es_bump2[1]]) - es_exact)
r1 = np.polyfit(np.log(sqrt_dof), np.log(err_es_p1), 1)[0]
plt.loglog(sqrt_dof, err_es_p1, marker='s', label="p=1, k={:.2f}".format(r1))

p = 2; Np = (p + 1) * (p + 2) / 2
sqrt_dof = np.sqrt(np.array([100 * Np, 400 * Np, 1600 * Np]))
err_es_p2 = np.abs(np.array([es_bump0[2], es_bump1[2], es_bump2[2]]) - es_exact)
r2 = np.polyfit(np.log(sqrt_dof[0:2]), np.log(err_es_p2[0:2]), 1)[0]
plt.loglog(sqrt_dof, err_es_p2, marker='s', label="p=2, k={:.2f}".format(r2))

plt.xlabel(r"$\sqrt{dof}$")
plt.ylabel(r"$Error\ of\ E_s$")
plt.legend(loc=3)

plt.savefig("conv.eps", bbox_inches='tight')

