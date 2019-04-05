#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(16, 4))
pdist0 = np.loadtxt("./bump2_p0/pressure_distribution.dat")

plt.subplot(1,3,1)
matplotlib.rcParams.update({'font.size': 10})
plt.scatter(pdist0[:, 0], -pdist0[:, 1], marker='.', s = 1)
plt.xlabel("x")
plt.ylabel(r"$-c_p$")
plt.title(r"$p = 0$")

pdist1 = np.loadtxt("./bump2_p1/pressure_distribution.dat")

plt.subplot(1,3,2)
plt.scatter(pdist1[:, 0], -pdist1[:, 1], marker='.', s = 1)
plt.xlabel("x")
plt.title(r"$p = 1$")


pdist2 = np.loadtxt("./bump2_p2/pressure_distribution.dat")

plt.subplot(1,3,3)
plt.scatter(pdist2[:, 0], -pdist2[:, 1], marker='.', s = 1)
plt.xlabel("x")
plt.title(r"$p = 2$")

plt.figure()
dp = pdist2 - pdist1
plt.scatter(pdist1[:, 0], dp[:, 1], marker='.', s = 1)


#plt.savefig("pdist.eps", bbox_inches="tight")

