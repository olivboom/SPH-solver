from sph_stub import SPH_main


def check_grid():
    domain = SPH_main()
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

def check_speed():
    domain = SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)

    for value in domain.particle_list:
        assert value.v[0] == 0
        assert value.v[1] == 0

def check_density():
    domain = SPH_main()
    domain.set_values()
    domain.initialise_grid()
    domain.place_points(domain.min_x, domain.max_x)

    for value in domain.particle_list:
        assert value.rho == 0
