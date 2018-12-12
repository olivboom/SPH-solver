import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import matplotlib
from sph_stub import SPH_main
from sph_stub import SPH_particle


def read_file_plot(solutions, option=None):
    '''Reads the output file from the numerical simulator. The x and y values
   (coordinates) of the particles from the simulator are stored inside two
   lists.
   Parameters:
       filename: the name of the file and extension to be read
       option: Default None: reads all particle parameters, including x and y
               coordinates, x and y velocities, pressure and density (rho)
               If option = 1: reads x and y coordinates and x velocity
               If option = 2: reads x and y coordinates and y velocity
               If option = 3: reads x and y coordinates and pressure
               If option = 4: reads x and y coordinates and density (rho)
   Returns: scatterplot of xcoordinate and ycoordinate and'''

    xcoord = []
    ycoord = []
    v_x = []
    v_y = []
    Pres = []
    rho = []
    output = []

    for i in range(len(solutions)):
        xcoord.append(solutions[i].x[0])
        ycoord.append(solutions[i].x[1])
        v_x.append(solutions[i].v[0])
        v_y.append(solutions[i].v[1])
        Pres.append(solutions[i].P)
        rho.append(solutions[i].rho)
    if option == 1:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        plt.scatter(xcoord, ycoord, c=v_x)
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-0.006, 0)
    if option == 2:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        plt.scatter(xcoord, ycoord, c=v_y)
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-100000 * 0.006, 0)
    if option == 3:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        plt.scatter(xcoord, ycoord, c=Pres)
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-0.006, 0)
    if option == 4:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
        plt.scatter(xcoord, ycoord, c=rho)
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-0.006, 0)
    if option is None:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 16))
        im1 = ax1.scatter(xcoord, ycoord, c=v_x)
        im2 = ax2.scatter(xcoord, ycoord, c=v_y)
        im3 = ax3.scatter(xcoord, ycoord, c=Pres)
        im4 = ax4.scatter(xcoord, ycoord, c=rho)
        fig.colorbar(im1, ax=ax1)
        fig.colorbar(im2, ax=ax2)
        fig.colorbar(im3, ax=ax3)
        fig.colorbar(im4, ax=ax4)

    plt.show()


def run():
    solutions = np.load('State.npy')

    for sol in solutions:
        read_file_plot(sol, option=2)

run()