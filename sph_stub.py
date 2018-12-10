"""SPH class to find nearest neighbours..."""

from itertools import count
import numpy as np


class SPH_main(object):
    """Primary SPH object"""

    def __init__(self):
        self.h = 0.0
        self.h_fac = 0.0
        self.dx = 0.0

        self.min_x = np.zeros(2)
        self.max_x = np.zeros(2)
        self.max_list = np.zeros(2, int)

        self.particle_list = []
        self.search_grid = np.empty((0, 0), object)

    def set_values(self):
        """Set simulation parameters."""

        self.min_x[:] = (0.0, 0.0)
        self.max_x[:] = (1.0, 1.0)
        self.dx = 0.02
        self.h_fac = 2
        self.h = self.dx * self.h_fac

    def initialise_grid(self):
        """Initalise simulation grid."""

        """Increases the minimum and maximum to account for the virtual particle padding that is required at boundaries"""
        self.min_x -= 2.0 * self.h
        self.max_x += 2.0 * self.h

        """Calculates the size of the array required to store the search array"""
        self.max_list = np.array((self.max_x - self.min_x) / (2.0 * self.h) + 1,
                                 int)

        self.search_grid = np.empty(self.max_list, object)

    def place_points(self, xmin, xmax):
        """Place points in a rectangle with a square spacing of size dx"""

        x = np.array(xmin)

        while x[0] <= xmax[0]:
            x[1] = xmin[1]
            while x[1] <= xmax[1]:
                particle = SPH_particle(self, x)
                particle.calc_index()
                if x[1] < 0.2:
                    self.particle_list.append(particle)
                elif x[1] < 0.5 and x[0] < 0.15:
                    self.particle_list.append(particle)
                elif x[1] > 0.5 and x[0] < 0:
                    self.particle_list.append(particle)
                elif x[1] > 0.2 and x[0] > 1:
                    self.particle_list.append(particle)
                elif x[1] > 1:
                    self.particle_list.append(particle)
                x[1] += self.dx
                # print(x[1])
            x[0] += self.dx

    def allocate_to_grid(self):
        """Allocate all the points to a grid in order to aid neighbour searching"""
        for i in range(self.max_list[0]):
            for j in range(self.max_list[1]):
                self.search_grid[i, j] = []

        for cnt in self.particle_list:
            self.search_grid[cnt.list_num[0], cnt.list_num[1]].append(cnt)

    def neighbour_iterate(self, part):
        """Find all the particles within 2h of the specified particle"""
        for i in range(max(0, part.list_num[0] - 1),
                       min(part.list_num[0] + 2, self.max_list[0])):
            for j in range(max(0, part.list_num[1] - 1),
                           min(part.list_num[1] + 2, self.max_list[1])):
                for other_part in self.search_grid[i, j]:
                    if part is not other_part:
                        dn = part.x - other_part.x
                        dist = np.sqrt(np.sum(dn ** 2))
                        if dist < 2.0 * self.h:
                            """This is only for demonstration - Your code will need to do all the particle to particle calculations at this point rather than simply displaying the vector to the neighbour"""
                            print("id:", other_part.id, "dn:", dn)


class SPH_particle(object):
    """Object containing all the properties for a single particle"""
    _ids = count(0)

    def __init__(self, main_data=None, x=np.zeros(2)):
        self.id = next(self._ids)
        self.main_data = main_data
        self.x = np.array(x)
        self.v = np.zeros(2)
        self.a = np.zeros(2)
        self.D = 0
        self.rho = 1.225 #0.0
        self.P = 0.0
        self.m = 0.0

    def calc_index(self):
        """Calculates the 2D integer index for the particle's location in the search grid"""
        self.list_num = np.array((self.x - self.main_data.min_x) /
                                 (2.0 * self.main_data.h), int)


def dw_dr(r, h):
    """
    Function to give the gradient of the smoothing kernel,
    includes compact support
    :param r: the positional vector between particle i and j
    :param h: the characteristic smoothing lenght
    :return:
    The gradient of the scaling factor associated with the smoothing kernel
    """
    constant = 10 / (7 * np.pi * h ** 2)
    r_mag = np.sqrt(r[0] ** 2 + r[1] ** 2)
    q = r_mag / h
    if 1 >= q >= 0:
        value = constant * (-3 * r_mag + 9 * r_mag ** 2 / 4)
    elif q >= 1 and q <= 2:
        value = constant * (-3 * (2 - r_mag / h) ** 2 / (4 * h))
    else:
        value = 0
    return value


def w(r, h):
    """
    Using a smoothing kernel returns the quantity A at location
    x
    :param r: the positional vector between particle i and j
    :param h: the characteristic smoothing length
    :return:
    The scaling factor associated with the smoothing kernel
    """
    constant = 10 / (7 * np.pi * h ** 2)
    rm = np.sqrt(r[0] ** 2 + r[1] ** 2)
    q = rm / h
    if q <= 1 and q >= 0:
        value = constant * (1 - 1.5 * q ** 2 + 0.75 * q ** 3)
    elif q >= 1 and q <= 2:
        value = constant * 0.25 * (2 - q) ** 3
    else:
        value = 0
    return value


def dv_dt(P1, P2, rho_1, rho_2, rij, vij, h, m2, mu):
    """
    Calculates the change in velocity for a time step
    for a
    :return:
    """
    mag_r = np.sqrt(sum(rij ** 2))
    eij = rij / mag_r
    dv = - m2 * (P1 / rho_1 ** 2 + P2 / rho_2 ** 2) * dw_dr(rij, h) * eij + (
                mu * m2 + (rho_1 ** -2 + rho_2 ** -2) * dw_dr(rij, h) * vij / mag_r)
    return dv


def drho_dt(rij, vij, h, m2):
    """
    Calculates the change in density for a time step
    :return:
    """
    mag_r = np.sqrt(sum(rij ** 2))
    eij = rij / mag_r
    drho = m2 * dw_dr(rij, h) * np.dot(vij, eij)
    return drho


def pressure(rho, rho0, c0, gamma):
    """Calculates the pressure for a given density"""
    B = (rho0 * c0 ** 2 / gamma) * rho0
    pres = B * ((rho / rho0) ** gamma - 1)
    return pres


def forward_euler(state, P1, P2, rho_1, rho_2, mu, c0, dt):
    """
    Because the equations are not coupled you can just time
    step through
    :return:
    """
    # Want to check that the right indices are taken in the right places
    vij = np.array((state[0], state[1]))
    dv = dv_dt(P1, P2, rho_1, rho_2, rij, vij, h, m2, mu)

    state[0] = state[0] + dt * dv[0]
    state[1] = state[1] + dt * dv[1]
    state[2] = state[2] + dt * drho_dt(rij, vij, h, m2)
    state[3] = pressure(rho, rho0, c0, gamma)
    return state


# """Create a single object of the main SPH type"""
domain = SPH_main()

"""Calls the function that sets the simulation parameters"""
domain.set_values()
"""Initialises the search grid"""
domain.initialise_grid()

"""Places particles in a grid over the entire domain - In your code you will need to place the fluid particles in only the appropriate locations"""
domain.place_points(domain.min_x, domain.max_x)

"""This is only for demonstration only - In your code these functions will need to be inside the simulation loop"""
"""This function needs to be called at each time step (or twice a time step if a second order time-stepping scheme is used)"""
domain.allocate_to_grid()
"""This example is only finding the neighbours for a single partle - this will need to be inside the simulation loop and will need to be called for every particle"""
# domain.neighbour_iterate(domain.particle_list[100])

rij = np.array([1, 1])
vij = np.array([2, 3])
mag_r = np.sqrt(sum(rij ** 2))
m2 = 5
gamma = 1
rho = 1.2
rho0 = 1.225
e_ij = rij / mag_r
h = 1.3 * 5
dt = 1
mu = 0.001
c0 = 20

i = 0
j = 1
P1 = domain.particle_list[i].P
P2 = domain.particle_list[j].P
rho_1 = domain.particle_list[i].rho
rho_2 = domain.particle_list[j].rho

state = [domain.particle_list[i].v[0], domain.particle_list[i].v[1], domain.particle_list[i].rho, domain.particle_list[i].P]
print(state)
state = forward_euler(state, P1, P2, rho_1, rho_2, mu, c0, dt)
print(state)


for i in range(len(domain.particle_list)):
    print(i)
    for j in range(1, len(domain.particle_list)):
        print(j)
        state = forward_euler(state, P1, P2, rho_1, rho_2, mu, c0, dt)



# fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))

# for i in range(0, len(domain.particle_list)):
#     ax1.plot(domain.particle_list[i].x[0], domain.particle_list[i].x[1], 'bo')
#
# plt.show()

# domain.neighbour_iterate(domain.particle_list[100])

# x_value = []
# y_value = []

# for value in domain.particle_list:
#     x_value.append(value.x[0])
#     y_value.append(value.x[1])
#     print(x_value)
#     print(y_value)
# ax3 = plt.subplot(111)
# ax3.plot(x_value, y_value, 'b.')
