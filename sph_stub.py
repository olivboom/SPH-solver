"""SPH class to find nearest neighbours..."""

from itertools import count
import matplotlib.pyplot as plt
import numpy as np
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
        self.normal = []

        self.min_x = np.zeros(2)
        self.max_x = np.zeros(2)
        self.max_list = np.zeros(2, int)

        self.particle_list = []  # np.array([])
        self.search_grid = np.empty((0, 0), object)
        self.log = []
        self.boundary_width = 0


    def set_values(self):
        """Set simulation parameters."""
        # self.scale = 30
        self.min_x[:] = (0.0, 0.0)
        self.max_x[:] = (5.0, 10.0)
        self.dx = 0.2 #0.02
        self.h_fac = 2
        self.h = self.dx * self.h_fac

        # Added quantities
        self.mu = 0.001
        self.rho0 = 1000
        self.g = np.array((0, -9.81))
        self.c0 = 20
        self.gamma = 7
        self.dt = 0.2# * self.h / self.c0
        self.B = (self.rho0 * self.c0 ** 2) / self.gamma
        self.boundary_width = 3

        self.t_cfl = 1
        self.t_f = 1
        self.t_A= 1
        self.dt_default = 0.1 * self.h / self.c0

        self.t_max = 30

    def initialise_grid(self):
        """Initalise simulation grid."""
        """Calculates the size of the array required to store the search array"""
        self.max_list = np.array((2 * self.boundary_width + self.max_x - self.min_x) / (2.0 * self.h) + 1,
                                 int)

        self.search_grid = np.empty(self.max_list, object)


    def place_points(self, xmin, xmax, geometry='default'):
        """Place points in a rectangle with a square spacing of size dx"""

        numx = int((xmax[0] - xmin[0]) / self.dx) - 1
        numy = int((xmax[1] - xmin[1]) / self.dx) - 1
        x_arr = np.linspace(xmin[0] + self.dx, xmax[0] - self.dx, numx)
        y_arr = np.linspace(xmin[1] + self.dx, xmax[1] - self.dx, numy)
        initial = np.empty((numx, numy))

        if geometry == 'default':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 2:
                        initial[i, j] = 2
                    elif y < 5 and x < 3:
                        initial[i, j] = 2

            # Anticlockwise from left wall
            self.normal = [np.array([1, 0]), np.array([0, 1]), np.array([-1, 0]), np.array([0, -1])]

        elif geometry == 'wave':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 8 and x < 3:
                        initial[i, j] = 2
                    elif y < 0.5 * x - 1.5 and x > 3 and x < 17:
                        initial[i, j] = 1
                    elif y < 7 and x > 17:
                        initial[i, j] = 1

        elif geometry == 'wave_2':
            for i, x in enumerate(x_arr):
                for j, y in enumerate(y_arr):
                    if y < 8 and x < 4:
                        initial[i, j] = 2
                    elif y < 0.5 * x - 4 and x > 8 and x < 16:
                        initial[i, j] = 1
                    elif y < 4 and x > 16:
                        initial[i, j] = 1

        initial = np.pad(initial, (self.boundary_width,
                                   self.boundary_width), mode='constant', constant_values=1)

        x_arr = np.linspace(xmin[0] - (self.boundary_width - 1) * self.dx,
                            xmax[0] - (self.boundary_width - 1) * self.dx, numx + 6)

        y_arr = np.linspace(xmin[1] - (self.boundary_width - 1) * self.dx,
                            xmax[1] - (self.boundary_width - 1) * self.dx, numy + 6)

        for i, x in enumerate(x_arr):
            for j, y in enumerate(y_arr):
                particle = SPH_particle(self, np.array([x, y]), m=self.dx ** 2 * self.rho0, P=0)
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


    def neighbour_iterate(self, part, wall_forcing='Leonard_Jones'):
        """Find all the particles within 2h of the specified particle
        and calculates the acceleration vector, density change and
        pressure at the new density
        """
        part.a = self.g
        part.D = 0.0
        if part.rho < 500:
            part.rho = 500
        if part.rho > 1500:
            part.rho = 1500

        part.can_see_wall = False
        for i in range(max(0, part.list_num[0] - 1),
                       min(part.list_num[0] + 2, self.max_list[0])):
            for j in range(max(0, part.list_num[1] - 1),
                           min(part.list_num[1] + 2, self.max_list[1])):
                for other_part in self.search_grid[i, j]:
                    if part is not other_part:
                        r_ij = part.x - other_part.x
                        mag_r_ij = np.sqrt((r_ij[0] ** 2) + (r_ij[1] ** 2))
                        if mag_r_ij < 2 * self.h:
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

                            if np.all(v_ij, 0):
                                t_cfl = self.h / np.sqrt(v_ij.dot(v_ij))

                                if t_cfl < self.t_cfl:
                                    self.t_cfl = t_cfl

                        if mag_r_ij < 2 * self.h and other_part.boundary is True and not part.boundary:
                            part.can_see_wall = True

        t_f = np.sqrt(self.h / np.sqrt(np.sum(part.a ** 2)))
        t_A = self.h / self.c0 * np.sqrt((part.rho / self.rho0) ** (self.gamma - 1))

        if t_f < self.t_f:
            self.t_f = t_f

        if t_A < self.t_A:
            self.t_A = t_A

        if wall_forcing == 'Leonard_Jones':

            if part.can_see_wall is True:
                for count, wall_normal in enumerate(self.normal):
                    if count < 2:
                        # Check sign on this
                        dist = part.x - self.min_x
                        perp_dist = np.dot(dist, wall_normal)
                        # print('distance:', perp_dist)
                    else:
                        dist = part.x - self.max_x
                        perp_dist = np.dot(dist, wall_normal)
                        # print('Location:', part.x)
                        # print('Distance', dist)
                        # print('Perp Dist', perp_dist)

                    d0 = 0.9 * self.dx
                    # print('q:', q)
                    q = perp_dist / d0
                    if q < 1:
                        if q < 0.1:
                            q = 0.1
                            part.v = np.array([0, 0])

                        fact = 1 / q
                        P_ref = (self.rho0 * self.c0 ** 2 / self.gamma) * ((1.05 ** 2) - 1)
                        factor = (fact ** 4 - fact ** 2) / perp_dist
                        acc_factor = wall_normal * factor * (P_ref / part.rho)

                        part.a = part.a + acc_factor

        elif wall_forcing == 'Artificial Forcing':

            if part.can_see_wall is True:
                for count, wall_normal in enumerate(self.normal):
                    if count < 2:
                        dist = part.x - self.min_x
                        perp_dist = np.dot(dist, wall_normal)
                    else:
                        dist = part.x - self.max_x
                        perp_dist = np.dot(dist, wall_normal)

                    d0 = 0.9 * self.dx
                    q = perp_dist / d0
                    if q < 1:
                        if q < 0.1:
                            q = 0.1
                        acc_factor = wall_normal * abs(part.a)

                        part.a = part.a + acc_factor
        else:
            pass

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
                        # print('q:', q)
                        # print('w:', w(q, self.h))
                        # print('Numerator', numerator)
                        if other_part.rho > 1500:
                            other_part.rho = 1500
                        denominator = denominator + (w(q, self.h) / other_part.rho)
        if denominator > 0:
            part.rho = numerator / denominator
        else:
            part.rho = 500
            print('Here')

        if part.rho < 500:
            part.rho = 500
        if part.rho > 1500:
            part.rho = 1500
            print('There')

    def forward_wrapper(self):
        """Stepping through using forwards Euler"""
        t = 0
        i = 0
        j = 0
        obj = []
        t_tracker = 0.05
        while t < self.t_max:
            t_tracker = t_tracker + self.dt
            if t_tracker > 0.05:
                with open('State.npy', 'wb') as fp:
                    pickle.dump(self.particle_list, fp)
                with open('State.npy', 'rb') as fp:
                    current = pickle.load(fp)
                obj.append(current)
                t_tracker = 0

            i = i + 1
            j = j + 1

            if i == 10:
                print('dt:', self.dt)
                ns.run([self.particle_list])
                i = 0

            if j == 10:
                print('Smoothing')
                for particle in self.particle_list:
                    self.density_smoothing(particle)
                j = 0

            t_in_1 = time.time()
            for particle in self.particle_list:
                self.neighbour_iterate(particle)

            # The 0.1 can be in the range of 0.1-0.3
            self.dt = 0.3 * min(self.t_cfl, self.t_f, self.t_A)
            if self.dt < self.dt_default:
                self.dt = self.dt_default

            t_out_1 = time.time()

            for particle in self.particle_list:
                particle.update_values(self.B, self.rho0, self.gamma, self.dt, self.min_x, self.max_x)

            # print('Neighbour Loop', (t_out_1 - t_in_1))
            # print('Time:', t)
            # print('Change in time:', self.dt)
            t += self.dt

        with open('State.npy', 'wb') as fp:
            pickle.dump(self.particle_list, fp)
        with open('State.npy', 'rb') as fp:
            current = pickle.load(fp)
        obj.append(current)
        with open('State.npy', 'wb') as fp:
            pickle.dump(obj, fp)

        # ns.run([self.particle_list])


class SPH_particle(object):
    """Object containing all the properties for a single particle"""
    _ids = count(0)

    def __init__(self, main_data=None, x=np.zeros(2), m=0, P=0):
        self.id = next(self._ids)
        self.main_data = main_data
        self.x = np.array(x)
        self.v = np.zeros(2)
        self.a = np.array([0, 0])
        self.D = 0
        self.rho = 1000
        self.P = P
        self.m = m
        self.boundary = False
        self.can_see_walls = False


    def calc_index(self):
        """Calculates the 2D integer index for the particle's location in the search grid"""
        self.list_num = np.array((self.x - self.main_data.min_x) /
                                 (2.0 * self.main_data.h), int)


    def update_values(self, B, rho0, gamma, dt, min_x, max_x, bounce=1):
        """Updates the state of the particle for one time step forwards"""
        if not self.boundary:
        #     new_x = self.x + (self.v * dt)
        #     if new_x[0] < min_x[0] or new_x[0] > max_x[0]:
        #         self.v = [bounce, 1] * self.v  #np.array([0, 0])
        #     elif new_x[1] < min_x[1] or new_x[1] > max_x[1]:
        #         self.v = [1, bounce] * self.v
        #         new_x = self.x + (self.v * dt)
        #     else:
        #         pass
            self.x = self.x + (self.v * dt)
            self.v = self.v + (self.a * dt)
                # self.x = new_x

        self.rho = self.rho + (self.D * dt)

        if self.rho < 500:
            self.rho = 500
            print('Set 500')
        elif self.rho > 1500:
            self.rho = 1500
            print('Set 1500')

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
