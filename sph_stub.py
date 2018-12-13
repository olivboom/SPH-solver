"""SPH class to find nearest neighbours..."""

from itertools import count
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from copy import deepcopy as copy
import time
import matplotlib.pylab as plty
from matplotlib import animation
import numpy_save as ns
import pickle


class SPH_main(object):
    """Primary SPH object"""

    def __init__(self):
        self.h = 0.0
        self.h_fac = 0.0
        self.dx = 0.0
        self.k = 0

        self.min_x = np.zeros(2)
        self.max_x = np.zeros(2)
        self.max_list = np.zeros(2, int)

        self.particle_list = []  # np.array([])
        self.search_grid = np.empty((0, 0), object)
        self.log = []


    def set_values(self):
        """Set simulation parameters."""
        # self.scale = 30
        self.min_x[:] = (0.0, 0.0)
        self.max_x[:] = (20.0, 10.0)
        self.dx = 0.5 #0.02
        self.h_fac = 2
        self.h = self.dx * self.h_fac
        self.k = 0.8

        # Added quantities
        self.mu = 0.001
        self.rho0 = 1000
        self.g = np.array((0, -9.81))  # 100000 factor to see difference in couple of small time steps
        self.c0 = 20
        self.gamma = 7
        self.dt = 0.1 * self.h / self.c0
        self.B = (self.rho0 * self.c0 ** 2) / self.gamma

        self.t_max = 200 * self.dt

    def initialise_grid(self):
        """Initalise simulation grid."""

        """Increases the minimum and maximum to account for the virtual particle padding that is required at boundaries"""
        self.min_x -= 2.0 * self.h
        self.max_x += 2.0 * self.h

        """Calculates the size of the array required to store the search array"""
        self.max_list = np.array((self.max_x - self.min_x) / (2.0 * self.h) + 1,
                                 int)

        self.search_grid = np.empty(self.max_list, object)


    def place_points(self, xmin, xmax, geometry='default'):
        """Place points in a rectangle with a square spacing of size dx"""

        x = np.array(xmin)
        numx = int((xmax[0] - xmin[0]) / self.dx)
        numy = int((xmax[1] - xmin[1]) / self.dx)
        x_arr = np.linspace(0, 20, numx)
        y_arr = np.linspace(0, 10, numy)
        initial = np.empty((numx, numy))

        if geometry == 'default':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 2:
                        initial[i, j] = 2
                    elif y < 5 and x < 3:
                        initial[i, j] = 2
            initial = np.pad(initial, (3, 3), mode='constant', constant_values=1)

        elif geometry == 'wave':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 8 and x < 3:
                        initial[i, j] = 2
                    elif y < 0.5 * x - 1.5 and x > 3 and x < 17:
                        initial[i, j] = 1
                    elif y < 7 and x > 17:
                        initial[i, j] = 1
            initial = np.pad(initial, (3, 3), mode='constant', constant_values=1)

        elif geometry == 'wave_2':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 8 and x < 4:
                        initial[i, j] = 2
                    elif y < 0.5 * x - 4 and x > 8 and x < 16:
                        initial[i, j] = 1
                    elif y < 4 and x > 16:
                        initial[i, j] = 1
            initial = np.pad(initial, (3, 3), mode='constant', constant_values=1)


        x_arr = np.linspace(xmin[0], xmax[0], numx + 6)
        y_arr = np.linspace(xmin[1], xmax[1], numy + 6)

        for i, x in enumerate(x_arr):
            for j, y in enumerate(y_arr):
                # TEMPORARY
                prefactor = (1000 / self.rho0) ** self.gamma
                P = (prefactor - 1) * self.B
                particle = SPH_particle(self, np.array([x, y]), m=self.dx ** 2 * self.rho0, P=P)
                particle.calc_index()
                if initial[i, j] == 2:
                    self.particle_list.append(particle)
                elif initial[i, j] == 1:
                    particle.boundary = True
                    self.particle_list.append(particle)


    def allocate_to_grid(self):
        """Allocate all the points to a grid in order to aid neighbour searching"""
        for i in range(self.max_list[0]):
            for j in range(self.max_list[1]):
                self.search_grid[i, j] = []

        for cnt in self.particle_list:
            self.search_grid[cnt.list_num[0], cnt.list_num[1]].append(cnt)


    def neighbour_iterate(self, part):
        """Find all the particles within 2h of the specified particle
        and calculates the acceleration vector, density change and
        pressure at the new density
        """
        # print(part.list_num)
        part.a = self.g
        part.D = 0.0
        for i in range(max(0, part.list_num[0] - 1),
                       min(part.list_num[0] + 2, self.max_list[0])):
            for j in range(max(0, part.list_num[1] - 1),
                           min(part.list_num[1] + 2, self.max_list[1])):
                for other_part in self.search_grid[i, j]:
                    if part is not other_part:
                        r_ij = part.x - other_part.x
                        mag_r_ij = np.sqrt((r_ij[0] ** 2) + (r_ij[1] ** 2))
                        if mag_r_ij < 2 * self.h:# and other_part.boundary == False:
                            mj = other_part.m
                            q = mag_r_ij / self.h
                            dwdr = dw_dr(q, self.h)
                            pre_fac = mj * dwdr / mag_r_ij

                            inv_rhoi2 = part.rho ** -2
                            inv_rhoj2 = other_part.rho ** -2
                            v_ij = part.v - other_part.v
                            fac2_1 = self.mu * (inv_rhoi2 + inv_rhoj2) * v_ij
                            Pi = part.P
                            Pj = other_part.P
                            fac2_2 = ((Pi * inv_rhoi2) + (Pj * inv_rhoj2)) * r_ij
                            post_fac = fac2_1 - fac2_2

                            part.a = part.a + (pre_fac * post_fac)
                            part.D = part.D + pre_fac * np.dot(v_ij, r_ij)

                        # elif mag_r_ij < self.k and other_part.boundary == True:
                        #
                        #     pre_factor = (self.c0 ** 2) * (self.k ** (self.gamma - 1) / (self.gamma * self.k)) / mag_r_ij
                        #     term_1 = (self.k / mag_r_ij) ** 2
                        #     acc_mag = pre_factor * (term_1 * (term_1 - 1))
                        #     part.a = part.a + (r_ij/mag_r_ij) * acc_mag

                            # if part.id == 200:
                            #     print("id:", part.id)
                            #     print('X coord:', part.x)
                            #     print('Velocity:', part.v)
                            #     print('Boundary:', part.boundary)
                            #     print('Acceleration', part.a)
                            #     print('Density change', part.D)


    def density_smoothing(self, part):
        numerator = 0
        denominator = 0
        for i in range(max(0, part.list_num[0] - 1),
                       min(part.list_num[0] + 2, self.max_list[0])):
            for j in range(max(0, part.list_num[1] - 1),
                           min(part.list_num[1] + 2, self.max_list[1])):
                for other_part in self.search_grid[i, j]:
                    r_ij = part.x - other_part.x
                    mag_r_ij = np.sqrt((r_ij[0] ** 2) + (r_ij[1] ** 2))
                    if mag_r_ij < 2 * self.h:
                        q = mag_r_ij / self.h
                        numerator = numerator + w(q, self.h)
                        denominator = denominator + (w(q, self.h) / other_part.rho)
        part.rho = numerator / denominator


    def forward_wrapper(self):
        """Stepping through using forwards Euler"""
        t = 0
        i = 0
        j = 0
        while t < self.t_max:
            i += 1
            j += 1
            # self.log.append(copy(self.particle_list))
            if i == 5:
                ns.run([self.particle_list])
                i = 0

            if j == 20:
                print('Smoothing')
                for particle in self.particle_list:
                    self.density_smoothing(particle)
                j = 0

            # plot the domain
            t_in_1 = time.time()
            for particle in self.particle_list:
                self.neighbour_iterate(particle)
            t_out_1 = time.time()

            print('Neighbour Loop', (t_out_1 - t_in_1))

            t_in = time.time()
            for particle in self.particle_list:
                particle.update_values(self.B, self.rho0, self.gamma, self.dt)
                # if particle.boundary is True:
                #     print('Boundary Coordinates: ', particle.x)

            t_out = time.time()

            t += self.dt

        ns.run([self.particle_list])
        # self.log.append(copy(self.particle_list))
        self.state_save()

    def state_save(self):
        np.save('State', self.log)


class SPH_particle(object):
    """Object containing all the properties for a single particle"""
    _ids = count(0)

    def __init__(self, main_data=None, x=np.zeros(2), m=0, P=0):
        self.id = next(self._ids)
        self.main_data = main_data
        self.x = np.array(x)
        self.v = np.zeros(2)
        self.a = np.array([0, -9.81])
        self.D = 0
        self.rho = 1000
        self.P = P
        self.m = m
        self.boundary = False


    def calc_index(self):
        """Calculates the 2D integer index for the particle's location in the search grid"""
        self.list_num = np.array((self.x - self.main_data.min_x) /
                                 (2.0 * self.main_data.h), int)

    def update_values(self, B, rho0, gamma, dt):
        """Updates the state of the particle for one time step forwards"""
        if self.boundary == False:
            self.x = self.x + (self.v * dt)
            self.v = self.v + (self.a * dt)

            # if (self.x[1]<self.main_data.h):
            #     print(self.a,' ',self.P,' ',self.rho)

        self.rho = self.rho + (self.D * dt)
        prefactor = self.rho / rho0
        self.P = (prefactor ** gamma - 1) * B


def dw_dr(q, h):
    """
    Function to give the gradient of the smoothing kernel,
    includes compact support
    :param r: the positional vector between particle i and j
    :param h: the characteristic smoothing lenght
    :return:
    The gradient of the scaling factor associated with the smoothing kernel
    """
    if q < 1:
        value = 10 / (7 * np.pi * h ** 3) * (-3 * q + 2.25 * q ** 2)
    else:
        value = -10 / (7 * np.pi * h ** 3) * 0.75 * (2 - q) ** 2
    return value


def w(q, h):
    """
    Using a smoothing kernel returns the quantity A at location
    x
    :param r: the positional vector between particle i and j
    :param h: the characteristic smoothing length
    :return:
    The scaling factor associated with the smoothing kernel
    """
    if q <= 1 and q >= 0:
        value = 10 / (7 * np.pi * h ** 2) * (1 - 1.5 * q ** 2 + 0.75 * q ** 3)
    elif q >= 1 and q <= 2:
        value = 10 / (7 * np.pi * h ** 2) * 0.25 * (2 - q) ** 3
    else:
        value = 0
    return value


if __name__ == "__main__":
    t0 = time.time()
    # """Create a single object of the main SPH type"""
    domain = SPH_main()

    domain.set_values()

    domain.initialise_grid()

    """Places particles in a grid over the entire domain"""
    domain.place_points(domain.min_x, domain.max_x)

    """This is only for demonstration only - In your code these functions will need to be inside the simulation loop"""
    """This function needs to be called at each time step (or twice a time step if a second order time-stepping scheme is used)"""
    domain.allocate_to_grid()

    """This example is only finding the neighbours for a single partle - this will need to be inside the simulation loop and will need to be called for every particle"""
    domain.forward_wrapper()

    # import numpy_save as ns

    t1 = time.time()

    print('Time to run:', t1 - t0)
