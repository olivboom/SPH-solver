import sph as sph

min_x = (0.0, 0.0) # lower domain boundary
max_x = (20.0, 10.0) # upper domain boundary
dx = 0.2 # change in x for initial condition
h_fac = 1.3 # h factor
b_width = 3 # number of boundary layers
t_max = 30 # maximum time to simulate up to 
geometry = 'default' # selects geometry // 'default' is a rectangle // 'wave' and 'wave_2' are beachheads
w_force = 'else' # wall forcing with Lennard-Jones // 'else' turns it off // 'Leonard_Jones' turns it on // 'Artificial Forcing' introduces corrective measures into L_J
q_fac = 0.1 # region within which inter-particle proximity correction takes place
bounce = -0.5 # behaviour of particles upon collision with boundary

sph.run(min_x=min_x, max_x=max_x, dx=dx, h_fac=h_fac,b_width=b_width, t_max=t_max,
        geometry=geometry, w_force=w_force, q_fac=q_fac, bounce=bounce)

option = 3 # colour bars /// 1: vx, 2: vy, 3: density, 4: pressure
import postprocessor as post


