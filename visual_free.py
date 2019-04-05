#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

res_p0 = np.loadtxt("./bump0_free/residual_p0.log")
res_p1 = np.loadtxt("./bump0_free/residual_p1.log")
res_p2 = np.loadtxt("./bump0_free/residual_p2.log")

plt.figure(figsize=(10,5))
matplotlib.rcParams.update({'font.size': 16})
plt.plot(res_p0[:, 0], res_p0[:, 1], label=r"$p = 0$")
plt.plot(res_p1[:, 0], res_p1[:, 1], label=r"$p = 1$")
plt.plot(res_p2[:, 0], res_p2[:, 1], label=r"$p = 2$")
plt.xlabel("Iteration Step")
plt.ylabel(r"$|R|_\infty$")
plt.legend()
plt.savefig("residual_free.eps", bbox_inches="tight")

