#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib as mpl
import numpy as np

def readdata(path):
    nodes_raw = np.loadtxt(path+"/nodes.dat")
    states_raw = np.loadtxt(path+"/states.dat")
    info = np.loadtxt(path+"/info.dat"); NE = int(info[0]); p = int(info[1])
    if p == 0:
        Np = 3
    else:
        Np = int((p + 1) * (p + 2) / 2)
    nodes = np.zeros([NE, Np, 2])
    states = np.zeros([NE, Np, 4])
    for i in range(NE):
        nodes[i, :, :] = nodes_raw[(i * Np) : (i * Np + Np), :]
        states[i, :, :] = states_raw[(i * Np) : (i * Np + Np), :]
    return nodes, states, NE, p

nodes, states, NE, p= readdata("bump2_p2")

fig = plt.figure(figsize=(12, 3))
cmap = 'jet'
M = np.sqrt(states[:, :, 1]**2 + states[:, :, 2]**2) / states[:, :, 0]
ax1 = fig.add_axes([0.0, 0.0, 0.9, 1.0])
for i in range(NE):
    print(i)
    x = np.asarray(nodes[i, :, 0]); y = np.asarray(nodes[i, :, 1]);
    z = np.asarray(M[i, :])
    if p == 2:
        triangles = [[0, 1, 3], [1, 4, 3], [1, 2, 4], [3, 4, 5]]
        tri_elem = [[0, 2, 5]]
    else:
        triangles = [[0, 1, 2]]
        tri_elem = [[0, 1, 2]]
    triang = mtri.Triangulation(x, y, triangles)
    triang_elem = mtri.Triangulation(x, y, tri_elem)
    ax1.tricontourf(triang, z, 5, vmin=np.min(M), vmax=np.max(M), cmap=cmap)
    #ax1.triplot(triang_elem, 'k-', linewidth=0.01)

plt.xlabel("x")
plt.ylabel("y")

ax2 = fig.add_axes([0.91, 0.0, 0.01, 1.0])
norm = mpl.colors.Normalize(vmin=np.min(M), vmax=np.max(M))
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('Mach Number')

plt.savefig("mach_bump2_p2.eps", bbox_inches='tight')
