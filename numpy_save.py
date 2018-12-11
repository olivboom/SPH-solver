import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from sph_stub import SPH_main
from sph_stub import SPH_particle


def read_file_plot(filename, option = None, time=-1):
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
    
    solutions = np.load(filename)
    xhist = []
    yhist = []
    vxhist = []
    vyhist = []
    Preshist = []
    rhohist = []
    
    for i in range(solutions.shape[0]):
        xcoord = []
        ycoord = []
        v_x = []
        v_y = []
        Pres = []
        rho = []

        current = solutions[i]

        for i in range(len(current)):
            xcoord.append(current[i].x[0])
            ycoord.append(current[i].x[1])
            v_x.append(current[i].v[0])
            v_y.append(current[i].v[1])
            Pres.append(current[i].P)
            rho.append(current[i].rho)

        xhist.append(xcoord)
        yhist.append(ycoord)
        vxhist.append(v_x)
        vyhist.append(v_y)
        Preshist.append(Pres)
        rhohist.append(rho)

    if option == 1:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        im1 = ax1.scatter(xhist[0], yhist[0], c=vyhist[0], vmin=-500, vmax=0)
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1], vmin=-500, vmax=0)
        plt.colorbar(im1, ax = ax1)
        plt.colorbar(im2, ax = ax2)
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax2.set_xlim(-0.1, 1.1)
        ax2.set_ylim(-0.1, 1.1)
    if option == 2:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        ax1.scatter(xhist[0], yhist[0], c=vyhist[0])
        ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1])
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-100000 * 0.006, 0)
    if option == 3:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        ax1.scatter(xhist[0], yhist[0], c=Preshist[0])
        ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1])
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-0.006, 0)
    if option == 4:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        ax1.scatter(xhist[0], yhist[0], c=rhohist[0])
        ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1])
        plt.colorbar()
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)
        plt.clim(-0.006, 0)
    if option is None:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 16))
        im1 = ax1.scatter(xhist[-1], yhist[-1], c=vxhist[-1])
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1])
        im3 = ax3.scatter(xhist[-1], yhist[-1], c=Preshist[-1])
        im4 = ax4.scatter(xhist[-1], yhist[-1], c=rhohist[-1])
        fig.colorbar(im1, ax = ax1)
        fig.colorbar(im2, ax = ax2)
        fig.colorbar(im3, ax = ax3)
        fig.colorbar(im4, ax = ax4)
        plt.show()
    
    fig = plt.figure(figsize=(8,8))
    ax1 = plt.subplot(111)
    scat = ax1.scatter([], [], c=[])
    
    def data_gen(i, xhist, yhist, vyhist, scat):
        # make sure initial condition is in plot
        scat = ax1.scatter(xhist[i], yhist[i], c=vyhist[i])
        return scat,
    
    anim = animation.FuncAnimation(fig, data_gen, frames=np.arange(0, solutions.shape[0]),
                                   fargs=(xhist, yhist, vyhist, scat), interval=100, blit=True, repeat=False)
    
    return anim

solutions = read_file_plot('State.npy', 1, time=-1)

#def main():
    

#if __name__ == '__main__':
    #main()