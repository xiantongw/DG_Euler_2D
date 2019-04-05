#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

res_p0 = np.loadtxt("./bump0_p0/residual.log")
res_p1 = np.loadtxt("./bump0_p1/residual.log")
res_p2 = np.loadtxt("./bump0_p2/residual.log")

plt.figure(figsize=(10,5))
matplotlib.rcParams.update({'font.size': 16})
plt.semilogy(res_p0[:, 0], res_p0[:, 1], label=r"$p = 0$")
plt.semilogy(res_p1[:, 0], res_p1[:, 1], label=r"$p = 1$")
plt.semilogy(res_p2[:, 0], res_p2[:, 1], label=r"$p = 2$")
plt.xlabel("Iteration Step")
plt.ylabel(r"$|R|_\infty$")
plt.legend()
plt.savefig("residual_full.eps", bbox_inches="tight")

