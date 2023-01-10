"""

This module can be considered the "front end" of the simulation.
Here users can choose input values to run the simulation with. The 'name'
parameter will become the name of the simulation under which the data will
be saved in the folder 'Simulation Data'.

"""

import simulation as sim
import initial_conditions as ic


def main():

    # Input parameters
    D = 2
    dens = 4
    T = 0.5
    num_tsteps = 40000
    tstep = 1E-3
    method = "Verlet"
    number_unitcells = 4**3
    pair_correlation_bin_size = 0.02
    name = 'Temp_05_Dens_4_2D'

    # Initializing system
    init_pos, N, box_len = ic.fcc_lattice(D, number_unitcells, dens)
    init_vel = ic.init_velocity(N, D, T)

    # Time evolution
    sim.simulate(N, number_unitcells, D, T, init_pos, init_vel, box_len, dens, num_tsteps, tstep, pair_correlation_bin_size, method, name)

    return


if __name__ == '__main__':
    main()