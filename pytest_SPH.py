import sph_stub

def test_grid():
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    
    x_value = []
    y_value = []

    for value in domain.particle_list:
       x_value.append(value.x[0])
       y_value.append(value.x[1])

    assert len(x_value) != 0
    assert len(y_value) != 0

def test_speed():
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    
    for value in domain.particle_list:
        assert value.v[0] == 0
        assert value.v[1] == 0

def test_density():
    domain = sph_stub.SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)
    domain.allocate_to_grid()
    
    for value in domain.particle_list:
        assert value.rho == 0
