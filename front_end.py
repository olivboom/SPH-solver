import sph as sph
import post_processor as post
import matplotlib.pyplot as plt

min_x = (0.0, 0.0)
max_x = (20.0, 10.0)
dx = 1
h_fac = 1.3
b_width = 3
t_max = 5
geometry = 'default'
w_force = 'else'
q_fac = 0.1
bounce = -0.5

sph.run(min_x=min_x, max_x=max_x, dx=dx, h_fac=h_fac,b_width=b_width, t_max=t_max,
        geometry=geometry, w_force=w_force, q_fac=q_fac, bounce=bounce)

option = 3 # colour bars /// 1: vx, 2: vy, 3: density, 4: pressure


post.run(option)

