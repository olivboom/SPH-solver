#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:38:07 2018
@author: Yujie-Z
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from sph_stub import SPH_main
from sph_stub import SPH_particle
import glob
from numpy import linalg as la


def last_5chars(x):
   return(x[-5:])


def read_file():
   filename = []
   for file in glob.glob("*.npy"):
       filename.append(file)
       filename1 = sorted(filename, key=last_5chars)

   xhist1 = [[] for m in range(len(filename1))]

   for i in range(len(filename1)):
       init_file = np.load(filename1[i])
       xhist = []

       for j in range(len(init_file)):
           xcoord = []

           current = init_file[j]

           for k in range(len(current)):
               xcoord.append(current[k].x[0])

           xhist.append(xcoord)

       xhist1[i] = xhist
   return xhist1
y = read_file()

def error_solver():
   error = []
   xhist1 = read_file()
   x_init = np.array(xhist1[-1])
   for i in range(len(xhist1)-1):
       x_array = np.array(xhist1[i])
       err = la.norm(x_array - x_init) / np.sqrt(len(x_array))
       error.append(err)
   return error

#x = [0.8, 1.0]
#error = error_solver()
#
#fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))
#ax1.loglog(x, error, 'b.')
