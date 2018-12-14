import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from sph_stub import SPH_main
from sph_stub import SPH_particle


def read_file_plot(filename, option=None, image=False):
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
    
    xboundhist = []
    yboundhist = []
    vxboundhist = []
    vyboundhist = []
    Presboundhist = []
    rhoboundhist = []

    vxmin = 0
    vxmax = 0
    vymin = 0
    vymax = 0
    Presmin = 0
    Presmax = 0
    rhomin = 0
    rhomax = 0

    print(len(solutions))

    for i in range(len(solutions)):
        xcoord = []
        ycoord = []
        v_x = []
        v_y = []
        Pres = []
        rho = []
        
        boundx = []
        boundy = []
        boundxv = []
        boundyv = []
        boundPres = []
        boundrho = []

        current = solutions[i]

        for i in range(len(current)):
            if current[i].boundary is True:
                boundx.append(current[i].x[0])
                boundy.append(current[i].x[1])
                boundxv.append(current[i].v[0])
                boundyv.append(current[i].v[1])
                boundPres.append(current[i].P)
                boundrho.append(current[i].rho)
               
            else:
                xcoord.append(current[i].x[0])
                ycoord.append(current[i].x[1])
                v_x.append(current[i].v[0])
                v_y.append(current[i].v[1])
                Pres.append(current[i].P)
                rho.append(current[i].rho)
                
            vxmin = np.amin([vxmin, current[i].v[0]])
            vxmax = np.amax([vxmax, current[i].v[0]])
            vymin = np.amin([vymin, current[i].v[1]])
            vymax = np.amax([vymax, current[i].v[1]])
            Presmin = np.amin([Presmin, current[i].P])
            Presmax = np.amax([Presmax, current[i].P])
            rhomin = np.amin([rhomin, current[i].rho])
            rhomax = np.amax([rhomin, current[i].rho])

        xhist.append(xcoord)
        yhist.append(ycoord)
        vxhist.append(v_x)
        vyhist.append(v_y)
        Preshist.append(Pres)
        rhohist.append(rho)
        
        xboundhist.append(boundx)
        yboundhist.append(boundy)
        vxboundhist.append(boundxv)
        vyboundhist.append(boundyv)
        Presboundhist.append(boundPres)
        rhoboundhist.append(boundrho)

    if option == 1:
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
        im1 = ax1.scatter(xhist[0], yhist[0], s=10, c=vxhist[0],
                          cmap='coolwarm', vmin=vxmin, vmax=vxmax)
        im2 = ax1.scatter(xboundhist[0], yboundhist[0], s=50, c=vxboundhist[0],
                          marker='s', cmap='coolwarm', vmin=vxmin, vmax=vxmax)
        im3 = ax2.scatter(xhist[-1], yhist[-1], s=10, c=vxhist[-1],
                          cmap='coolwarm', vmin=vxmin, vmax=vxmax)
        im4 = ax2.scatter(xboundhist[-1], yboundhist[-1], s=50, c=vxboundhist[-1],
                          marker='s', cmap='coolwarm', vmin=vxmin, vmax=vxmax)
#        plt.colorbar(im1, ax=ax1)
#        plt.colorbar(im2, ax=ax2)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax2.set_xlim(-5, 25)
        ax2.set_ylim(-5, 15)

        if image is True:
            for i in range(len(solutions)):
                fig, ax = plt.subplots(1, 1, figsize=(16, 16))
                im = ax.scatter(xhist[i], yhist[i], s=10, c=vxhist[i],
                                cmap='coolwarm', vmin=vxmin, vmax=vxmax)
                im = ax.scatter(xboundhist[i], yboundhist[i],
                                s=50, c=vxboundhist[i], marker='s', cmap='coolwarm',
                                vmin=vxmin, vmax=vxmax)
#                plt.colorbar(im, ax=ax)
                ax.set_xlim(-5, 25)
                ax.set_ylim(-5, 15)
                plt.savefig('.\plot_images\image{}.png'.format(i))
                plt.close()

        fig3 = plt.figure(figsize=(16, 16))
        scat = ax1.scatter([], [], c=[], s=10, cmap='coolwarm',
                           vmin=vxmin, vmax=vxmax)
        scat2 = ax1.scatter([], [], c=[], s=50, marker='s', cmap='coolwarm',
                           vmin=vxmin, vmax=vxmax)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax1 = plt.subplot(111)

        def animate(i, xhist, yhist, xboundhist, yboundhist, vxhist, vxboundhist,
                    scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=vxhist[i], s=10,
                               cmap='coolwarm', vmin=vxmin, vmax=vxmax)
            scat2 = ax1.scatter(xboundhist[i], yboundhist[i],
                                c=vxboundhist[i], s=50,
                                marker='s', cmap='coolwarm', vmin=vxmin,
                                vmax=vxmax)
            ax1.set_xlim(-5, 25)
            ax1.set_ylim(-5, 15)
            return scat, scat2

        anim = animation.FuncAnimation(fig3, animate,
                                       frames=np.arange(0, len(solutions)),
                                       fargs=(xhist, yhist, xboundhist,
                                              yboundhist, vxhist, 
                                              vxboundhist, scat),
                                       interval=60, blit=True, repeat=True)

        return anim

    elif option == 2:
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
        im1 = ax1.scatter(xhist[0], yhist[0], s=10, c=vyhist[0],
                          cmap='coolwarm', vmin=vymin, vmax=vymax)
        im2 = ax1.scatter(xboundhist[0], yboundhist[0], s=50, c=vyboundhist[0],
                          marker='s', cmap='coolwarm', vmin=vymin, vmax=vymax)
        im3 = ax2.scatter(xhist[-1], yhist[-1], s=10, c=vyhist[-1],
                          cmap='coolwarm', vmin=vymin, vmax=vymax)
        im4 = ax2.scatter(xboundhist[-1], yboundhist[-1], s=50, c=vyboundhist[-1],
                          marker='s', cmap='coolwarm', vmin=vymin, vmax=vymax)
#        plt.colorbar(im1, ax=ax1)
#        plt.colorbar(im2, ax=ax2)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax2.set_xlim(-5, 25)
        ax2.set_ylim(-5, 15)

        if image is True:
            for i in range(len(solutions)):
                fig, ax = plt.subplots(1, 1, figsize=(16, 16))
                im = ax.scatter(xhist[i], yhist[i], s=10, c=vyhist[i],
                                cmap='coolwarm', vmin=vxmin, vmax=vxmax)
#                plt.colorbar(im, ax=ax)
                ax.set_xlim(-5, 25)
                ax.set_ylim(-5, 17)
                plt.savefig('.\plot_images\image{}.png'.format(i))
                plt.close()

        fig3 = plt.figure(figsize=(16, 16))
        scat = ax1.scatter([], [], c=[], s=10, cmap='coolwarm',
                           vmin=vymin, vmax=vymax)
        scat2 = ax1.scatter([], [], c=[], s=50, marker='s', cmap='coolwarm',
                           vmin=vymin, vmax=vymax)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax1 = plt.subplot(111)

        def animate(i, xhist, yhist, xboundhist, yboundhist, vyhist, vyboundhist,
                    scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=vyhist[i], s=10,
                               cmap='coolwarm', vmin=vymin, vmax=vymax)
            scat2 = ax1.scatter(xboundhist[i], yboundhist[i],
                                c=vyboundhist[i], s=50,
                                marker='s', cmap='coolwarm', vmin=vymin,
                                vmax=vymax)
            ax1.set_xlim(-5, 25)
            ax1.set_ylim(-5, 15)
            return scat, scat2

        anim = animation.FuncAnimation(fig3, animate,
                                       frames=np.arange(0, len(solutions)),
                                       fargs=(xhist, yhist, xboundhist,
                                              yboundhist, vyhist, 
                                              vyboundhist, scat),
                                       interval=60, blit=True, repeat=True)

        return anim

    elif option == 3:
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
        im1 = ax1.scatter(xhist[0], yhist[0], s=10, c=Preshist[0],
                          cmap='coolwarm', vmin=Presmin, vmax=Presmax)
        im2 = ax1.scatter(xboundhist[0], yboundhist[0], s=50, c=Presboundhist[0],
                          marker='s', cmap='coolwarm', vmin=Presmin, vmax=Presmax)
        im3 = ax2.scatter(xhist[-1], yhist[-1], s=10, c=Preshist[-1],
                          cmap='coolwarm', vmin=Presmin, vmax=Presmax)
        im4 = ax2.scatter(xboundhist[-1], yboundhist[-1], s=50, c=Presboundhist[-1],
                          marker='s', cmap='coolwarm', vmin=Presmin, vmax=Presmax)
#        plt.colorbar(im1, ax=ax1)
#        plt.colorbar(im2, ax=ax2)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax2.set_xlim(-5, 25)
        ax2.set_ylim(-5, 15)

        if image is True:
            for i in range(len(solutions)):
                fig, ax = plt.subplots(1, 1, figsize=(16, 16))
                im = ax.scatter(xhist[i], yhist[i], s=10, c=Preshist[i],
                                cmap='coolwarm', vmin=Presmin, vmax=Presmax)
#                plt.colorbar(im, ax=ax)
                ax.set_xlim(-5, 25)
                ax.set_ylim(-5, 17)
                plt.savefig('.\plot_images\image{}.png'.format(i))
                plt.close()

        fig3 = plt.figure(figsize=(16, 16))
        scat = ax1.scatter([], [], c=[], s=10, cmap='coolwarm',
                           vmin=Presmin, vmax=Presmax)
        scat2 = ax1.scatter([], [], c=[], s=50, marker='s', cmap='coolwarm',
                           vmin=Presmin, vmax=Presmax)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax1 = plt.subplot(111)

        def animate(i, xhist, yhist, xboundhist, yboundhist, Preshist, Presboundhist,
                    scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=Preshist[i], s=10,
                               cmap='coolwarm', vmin=Presmin, vmax=Presmax)
            scat2 = ax1.scatter(xboundhist[i], yboundhist[i],
                                c=Presboundhist[i], s=50,
                                marker='s', cmap='coolwarm', vmin=Presmin,
                                vmax=Presmax)
            ax1.set_xlim(-5, 25)
            ax1.set_ylim(-5, 15)
            return scat, scat2

        anim = animation.FuncAnimation(fig3, animate,
                                       frames=np.arange(0, len(solutions)),
                                       fargs=(xhist, yhist, xboundhist,
                                              yboundhist, Preshist, 
                                              Presboundhist, scat),
                                       interval=60, blit=True, repeat=True)

        return anim

    elif option == 4:
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
        im1 = ax1.scatter(xhist[0], yhist[0], s=10, c=rhohist[0],
                          cmap='coolwarm', vmin=rhomin, vmax=rhomax)
        im2 = ax1.scatter(xboundhist[0], yboundhist[0], s=50, c=rhoboundhist[0],
                          marker='s', cmap='coolwarm', vmin=rhomin, vmax=rhomax)
        im3 = ax2.scatter(xhist[-1], yhist[-1], s=10, c=rhohist[-1],
                          cmap='coolwarm', vmin=rhomin, vmax=rhomax)
        im4 = ax2.scatter(xboundhist[-1], yboundhist[-1], s=50, c=rhoboundhist[-1],
                          marker='s', cmap='coolwarm', vmin=rhomin, vmax=rhomax)
#        plt.colorbar(im1, ax=ax1)
#        plt.colorbar(im2, ax=ax2)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax2.set_xlim(-5, 25)
        ax2.set_ylim(-5, 15)

        if image is True:
            for i in range(len(solutions)):
                fig, ax = plt.subplots(1, 1, figsize=(16, 16))
                im = ax.scatter(xhist[i], yhist[i], s=10, c=rhohist[i],
                                cmap='coolwarm', vmin=rhomin, vmax=rhomax)
#                plt.colorbar(im, ax=ax)
                ax.set_xlim(-5, 25)
                ax.set_ylim(-5, 17)
                plt.savefig('.\plot_images\image{}.png'.format(i))
                plt.close()

        fig3 = plt.figure(figsize=(16, 16))
        scat = ax1.scatter([], [], c=[], s=10, cmap='coolwarm',
                           vmin=rhomin, vmax=rhomax)
        scat2 = ax1.scatter([], [], c=[], s=50, marker='s', cmap='coolwarm',
                           vmin=rhomin, vmax=rhomax)
        ax1.set_xlim(-5, 25)
        ax1.set_ylim(-5, 15)
        ax1 = plt.subplot(111)

        def animate(i, xhist, yhist, xboundhist, yboundhist, rhohist, rhoboundhist,
                    scat):
            # make sure initial condition is in plot
            scat = ax1.scatter(xhist[i], yhist[i], c=rhohist[i], s=10,
                               cmap='coolwarm', vmin=rhomin, vmax=rhomax)
            scat2 = ax1.scatter(xboundhist[i], yboundhist[i],
                                c=rhoboundhist[i], s=50,
                                marker='s', cmap='coolwarm', vmin=rhomin,
                                vmax=rhomax)
            ax1.set_xlim(-5, 25)
            ax1.set_ylim(-5, 15)
            return scat, scat2

        anim = animation.FuncAnimation(fig3, animate,
                                       frames=np.arange(0, len(solutions)),
                                       fargs=(xhist, yhist, xboundhist,
                                              yboundhist, rhohist, 
                                              rhoboundhist, scat),
                                       interval=60, blit=True, repeat=True)

        return anim


solutions = read_file_plot('State.npy', 3)

