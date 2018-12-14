"""
Created on Thu Dec 13 15:24:01 2018

@author: Yujie-Z
"""
import matplotlib.pyplot as plt
import numpy as np
from sph_stub import SPH_main
from sph_stub import SPH_particle
from numpy import linalg as la
import pandas as pd

def read_h(filename, dx):
    c0 = 20
    h_fac = 2
    h = dx * h_fac
    dt = 0.1 * h * c0
    solution = np.load(filename)
    xlist = []
    ylist = []
    for j in range(len(solution)):
        current = solution[j]
        h = []
        for i in range(len(current)):
            h1 = current[i].x
            if current[i].boundary == False and 0 < h1[0] < 20:
                h.append(h1)
        harray = np.array(h)
        L1 = pd.DataFrame(harray).groupby(0, as_index=False)[1].max().values.tolist()
        L = np.array(L1)
        ymax = np.max(L[:,1])
        index = np.where(L[:,1]==ymax)[0]
        max_value = L[index[int(len(index)/2)]]
        xlist.append(max_value[0])
        ylist.append(max_value[1])
    xmax = np.array(xlist)
    ymax = np.array(ylist)
    v_max = (xmax[1:] - xmax[0:-1])/dt
    v_max = np.abs(v_max)
    v_average = np.sum(v_max)/np.size(v_max)
    return v_average, dt, ymax

def analytical_solution(dx):
    g = 9.81
    h0 = (((20/dx)+1)*(3/20)*0.5 + ((20/dx)+1)*(17/20)*0.2)/((20/dx)+1)
    ana_solution = np.sqrt(g*h0)
    return ana_solution


dx = 0.5
y_ana = analytical_solution(dx)
y, dt, ymax = read_h('State.npy', dx)
err = np.abs(y - y_ana)
print('dt= ', dt)
print('dx= ', dx)
print('ymax', ymax)
print('discrepancy',err)
