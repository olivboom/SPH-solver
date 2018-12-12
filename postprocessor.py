import matplotlib.pyplot as plt
import matplotlib
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

    vxmin = 0
    vxmax = 0
    vymin = 0
    vymax = 0
    Presmin = 0
    Presmax = 0
    rhomin = 0
    rhomax = 0

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
            vxmin = np.amin([vxmin, current[i].v[0]])
            vxmax = np.amax([vxmax, current[i].v[0]])
            vymin = np.amin([vymin, current[i].v[0]])
            vymax = np.amax([vymax, current[i].v[0]])

            Pres.append(current[i].P)
            Presmin = np.amin([Presmin, current[i].P])
            Presmax = np.amax([Presmax, current[i].P])

            rho.append(current[i].rho)
            rhomin = np.amin([rhomin, current[i].rho])
            rhomax = np.amax([rhomin, current[i].rho])

        xhist.append(xcoord)
        yhist.append(ycoord)
        vxhist.append(v_x)
        vyhist.append(v_y)
        Preshist.append(Pres)
        rhohist.append(rho)

    if option == 1:
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        im1 = ax1.scatter(xhist[0], yhist[0], c=vxhist[0],
                          vmin=vxmin, vmax=vxmax)
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=vxhist[-1],
                          vmin=vxmin, vmax=vxmax)
        plt.colorbar(im1, ax = ax1)
        plt.colorbar(im2, ax = ax2)
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax2.set_xlim(-0.1, 1.1)
        ax2.set_ylim(-0.1, 1.1)
        

#        for i in range(solutions.shape[0]):
#            fig, ax = plt.subplots(1, 1, figsize=(8,8))
#            im = ax.scatter(xhist[i], yhist[i], c=vxhist[i],
#                                     vmin=vxmin, vmax=vxmax)
#            plt.colorbar(im, ax=ax)
#            plt.savefig('image{}.png'.format(i))
#            plt.close()
        
        fig3 = plt.figure(figsize=(8,8))
        scat = ax1.scatter(xhist[0], yhist[0],c=vxhist[0],
                           vmin=vxmin, vmax=vxmax)
        ax1 = plt.subplot(111)
        
#        def init():
#            scat = ax1.scatter(xhist[i], yhist[i], c=vxhist[i], vmin=vxmin,
#                               vmax=vxmax)
#            return scat,
        
        def animate(i, xhist, yhist, vxhist, scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=vxhist[i],
                                     vmin=vxmin, vmax=vxmax)
            return scat,

        anim = animation.FuncAnimation(fig3, animate,
                                       frames=np.arange(0, solutions.shape[0]),
                                       fargs=(xhist, yhist, vxhist, scat),
                                       interval=1, blit=True, repeat=True)
        
        return anim
    
    elif option == 2:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        im1 = ax1.scatter(xhist[0], yhist[0], c=vyhist[0],
                          vmin=vymin, vmax=vymax)
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=vyhist[-1],
                          vmin=vymin, vmax=vymax)
        plt.colorbar(im1, ax = ax1)
        plt.colorbar(im2, ax = ax2)
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax2.set_xlim(-0.1, 1.1)
        ax2.set_ylim(-0.1, 1.1)
        
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
    
    elif option == 3:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        im1 = ax1.scatter(xhist[0], yhist[0], c=Preshist[0], vmin=Presmin, vmax=Presmax)
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=Preshist[-1], vmin=Presmin, vmax=Presmax)
        plt.colorbar(im1, ax = ax1)
        plt.colorbar(im2, ax = ax2)
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax2.set_xlim(-0.1, 1.1)
        ax2.set_ylim(-0.1, 1.1)
        
        fig = plt.figure(figsize=(8,8))
        ax1 = plt.subplot(111)
        scat = ax1.scatter([], [], c=[])
        
        def data_gen(i, xhist, yhist, Preshist, scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=Preshist[i])
            return scat,
        
        anim = animation.FuncAnimation(fig, data_gen, frames=np.arange(0, solutions.shape[0]),
                                       fargs=(xhist, yhist, Preshist, scat), interval=100, blit=True, repeat=False)
        
        return anim
    
    elif option == 4:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        im1 = ax1.scatter(xhist[0], yhist[0], c=rhohist[0], vmin=rhomin, vmax=rhomax)
        im2 = ax2.scatter(xhist[-1], yhist[-1], c=rhohist[-1], vmin=rhomin, vmax=rhomax)
        plt.colorbar(im1, ax = ax1)
        plt.colorbar(im2, ax = ax2)
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax2.set_xlim(-0.1, 1.1)
        ax2.set_ylim(-0.1, 1.1)
        
        fig = plt.figure(figsize=(8,8))
        ax1 = plt.subplot(111)
        scat = ax1.scatter([], [], c=[])
        
        def data_gen(i, xhist, yhist, rhohist, scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=rhohist[i])
            return scat,
        
        anim = animation.FuncAnimation(fig, data_gen, frames=np.arange(0, solutions.shape[0]),
                                       fargs=(xhist, yhist, rhohist, scat), interval=100, blit=True, repeat=False)
        
        return anim

    

solutions = read_file_plot('State.npy', 1, time=-1)
