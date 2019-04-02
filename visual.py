#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 23:46:34 2019

@author: xtwang
"""

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

## Create triangulation.
#x = np.asarray([0, 1, 2, 3, 0.5, 1.5, 2.5, 1, 2, 1.5])
#y = np.asarray([0, 0, 0, 0, 1.0, 1.0, 1.0, 2, 2, 3.0])
#triangles = [[0, 1, 4], [1, 2, 5], [2, 3, 6], [1, 5, 4], [2, 6, 5], [4, 5, 7],
#             [5, 6, 8], [5, 8, 7], [7, 8, 9]]
#triang = mtri.Triangulation(x, y, triangles)
#
## Interpolate to regularly-spaced quad grid.
#z = np.cos(1.5 * x) * np.cos(1.5 * y)
#
## Set up the figure
#fig, axs = plt.subplots(nrows=1, ncols=1)
#
## Plot the triangulation.
#axs.tricontourf(triang, z, 20)
#axs.triplot(triang, 'ko-')
#axs.set_title('Triangular grid')
#
#fig.tight_layout()
#plt.show()


#def read_data():
#    nodes_raw = np.load("nodes.dat")
#    states_raw = np.load("states.dat")

def readdata():
    nodes_raw = np.loadtxt("nodes.dat")
    states_raw = np.loadtxt("states.dat")
    info = np.loadtxt("info.dat"); NE = int(info[0]); p = int(info[1])
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

nodes, states, NE, p= readdata()

plt.figure(figsize=(10, 3))
for i in range(NE):
    x = np.asarray(nodes[i, :, 0]); y = np.asarray(nodes[i, :, 1]);
    z = np.asarray(states[i, :, 0])
    if p == 2:
        triangles = [[0, 1, 3], [1, 4, 3], [1, 2, 4], [3, 4, 5]]
        tri_elem = [[0, 2, 5]]
    else:
        triangles = [[0, 1, 2]]
        tri_elem = [[0, 1, 2]]
    triang = mtri.Triangulation(x, y, triangles)
    triang_elem = mtri.Triangulation(x, y, tri_elem)
    
    plt.tricontourf(triang, z, 20, vmin=np.min(states[:, :, 0]), vmax=np.max(states[:, :, 0]))
    plt.triplot(triang_elem, 'k-', linewidth=0.2)