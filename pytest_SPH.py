import sph_stub

def test_grid():
    ''' Tests if the initial grid has been populated by particles
    '''
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    
    x_value = []
    y_value = []

    for particle in domain.particle_list:
       x_value.append(particle.x[0])
       y_value.append(particle.x[1])

    assert len(x_value) != 0
    assert len(y_value) != 0

def test_speed():
    ''' Tests if any particles have exceeded the speed of sound
    '''
    
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    domain.t_max = 2 * domain.dt
    domain.forward_wrapper()
    
    for time_log in domain.log:
        for particle in time_log:
            #assert particle.v[0] == 0
            #assert particle.v[1] == 0
            assert particle.v[0] <= domain.c0
            assert particle.v[1] <= domain.c0
    
    
def test_density():
    ''' Tests if a negative density has been given and whether any
    particle density has exceeded the reference density by a factor of 1.5
    '''
    
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    domain.t_max = 2 * domain.dt
    domain.forward_wrapper()
    
    for time_log in domain.log:
        for particle in time_log:    
            assert particle.rho > 0
            assert particle.rho<= 1.5 * domain.rho0

def test_mass_conserve():
    ''' Tests whether any particles have fallen out of the initial domain
    '''
    
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    x_min = domain.min_x[0]
    x_max = domain.max_x[0]
    y_min = domain.min_x[1]
    y_max = domain.min_x[1]
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    domain.t_max = 2 * domain.dt
    domain.forward_wrapper()
    
    for time_log in domain.log:
        for particle in time_log:
            assert x_min <= particle.x[0] <= x_max
            assert y_min <= particle.x[1] <= y_max

    