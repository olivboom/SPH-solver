# Team Atlantic - A Smoothed Particle Hydrodynamics Solver



## Feature

This project simulates a dam break scenario and solves the resultant fluid flow based on the Smoothed Particle Hydrodynamics (SPH) method.

The key outputs from the simulation program are:

-  The particle data (including the locations, velocity, acceleration, density, pressure, mass of each individual fluid and boundary particle) stored in  an output file after a predetermined end time.
-  A plot of the particle data produced by the SPH simulator at a specific time.
-  An animation of the simulation results.



# Preview

![](https://raw.github.com/meolu/walle-web/master/screenshot/projects.png)

...................need a figure to replace it

# Installation

To run the numerical simulator, download the module sph_stub.py。 Calling the module in the command prompt should automatically produce a ‘State.npy’ file into the current working directory of the user which can then be read by the post-processor.
Below is an example of how to run the simulator:

**`$ python sph_stub.py`**

* Note, in the example above, the simulation will be based on pre-existing parameters not given by the user.
* To create an animation of the results, download the module postprocessor.py . Ensure that sph_stub.py is already downloaded before trying to operate postprocessor.py.
* Within a python kernel, the user can create the animation by the following steps:
  - import postprocessor
  - create a variable
  - execute the following command: variable = postprocessor.read_file_plot(‘State.npy’, option, time=-1, image=False)
    - if option=1 : creates an animation of the location of particles and the x-direction velocity from the initial time to the final simulation time
    - if option=2  : creates an animation of the location of particles and the y-direction velocity from the initial time to the final simulation time
    - if option=3  : creates an animation of the location of particles and the Pressure from the initial time to the final simulation time
    - if option=4: creates an animation of the location of particles and the density from the initial time to the final simulation time
    - if image=True: if the user selects this parameter to be true, then the postprocessor will save a series of snapshots at each time step of the simulation and save them into the user’s current working directory



# Structure

## The sph_stub module:

The sph_stub module is mainly used for generating a  domain and dynamic particles (and their properties) in the meshed domain.

 ### The  SPH_main class: 
The SPH_main class, within "sph_stub.py", used to initialise a single object of the main SPH type (a initial domain containing boundaries and liquid particles) <br>

#### Class methods:

##### **SPH_main():**

- This is the class object that creates the domain onto which the particles will be placed

**set_values():**

- This method sets the simulation domain parameters.
- Domain reference values such as the fluid speed of sound, the domain reference density and the solution end time are given .<br>

**initialise_grid():**

- Initialises a simulation grid and calculates the size of the array needed for an attribute called search_grid[]. search_grid[] is used later on to calculate the contributions of neighbouring particles in neighbour_iterate().

**place_points( xmin,  xmax):**

- within the initialised grid, a particle at every increment of dx is created using the SPH_particle() class. The particles are assigned an index relative to the search grid and then appended into a list called particle_list[].

**allocate_to_grid()**

- particles are assigned to a search grid using the index created earlier.

**neighbour_iterate(part)**

* this function iterates through every particle in the domain. It quantifies the contributions of the neighbours by using it's search_grid index and then approximates the derivatives in the navier-stokes equation using a smoothing kernel

**forward_wrapper( )**

-  this function time steps the particles in the particle list. At every time increment, neighbour_iterate and update_values are used to update the particle parameters for the next time increment
-  at the end of the forward_wrapper(), the each particle_list at prespecified time increment is combined and then saved into a file called 'State.npy'.

### The SPH_particle Class:

The SPH_particle class, within "sph_stub.py" is a class which creates an object of a single particle with attributes. The class contains the attributes including position (x), velocity (v), acceleration (a), density (rho), pressure (P) and the instantaneous change in density (D).

Within the SPH_particle class are class functions:

**calc_index()**

- this function gives assigns a list number to the current particle. This list number is used to allocate the particle to the search grid

**update_values()**

- this function updates the SPH_particle class attributes using derivatives calculated in neighbour_iterate()

# Testing

Testing is performed using pytest.  Four sanity checks are implemented:

-  test_grid() : a domain is created using and particles are allocated to the grid. This test checks whether or not if particles have actually been stored inside the particle_list (whether or not the particles exist)
- test_speed() : a solutions file created from the simulator is loaded. At every time interval, for every particle list, the velocity of particle is checked to ensure that it does not exceed the domain speed of sound
- test_density() : a solutions file created from the simulator is loaded. At every time interval, for every particle list, the particle density is checked to ensure it is not negative and it does not exceed the reference density by a factor of 1.5
- test_mass_conserve() : a solutions file created from the simulator is loaded. At every time interval, for every particle list, the particle coordinate is checked to ensure that it is within the bounds of the minimum and maximum of the domain

# Roadmap




  # Contributors

  - [**Professor Stephen Neethling**](https://www.imperial.ac.uk/people/s.neethling)
  - [**Atlantic group**](https://github.com/msc-acse/acse-4-project-2-atlantic)
