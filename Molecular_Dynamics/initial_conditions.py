"""

This module contains the function initializing the systems by setting the
initial particle positions according to a center faced chrystal lattice
and velocities.

"""

import numpy as np
import math as m


def fcc_lattice(D, number_unitcells, dens):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    D : int
        The number of dimensions in the system
    number_unitcells : int
        Sets the number of unit cells
    dens : float
        The density of the system defined as N/V

    Returns
    -------
    init_pos : np.ndarray
        Array of particle coordinates
    N : int
        The number of particles in the system
    """

    init_pos = np.empty((0, D), float)
    if D == 2:
        # Determining number of steps in the x and y direction
        it = int(number_unitcells**(1/2))
        # Calculating the box length
        a = m.pow(2/dens, 1/2)
        box_len = a * it
        # Number of particles inside the box
        N = 0
        for i in range(it):
            for j in range(it):
                init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a, (j + 1/3)*a])))
                N += 1
                init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a + a/2, (j + 1/3)*a + a/2])))
                N += 1

    if D == 3:
        # Determining number of steps in the x and y direction
        it = int(round(number_unitcells**(1/3), 0))
        # Calculating the box length
        a = m.pow(4/dens, 1/3)
        box_len = a * it
        # Number of particles inside the box
        N = 0
        for i in range(it):
            for j in range(it):
                for k in range(it):
                    init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a, (j + 1/3)*a, (k + 1/3)*a])))
                    N += 1
                    init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a + a/2, (j + 1/3)*a + a/2, (k + 1/3)*a])))
                    N += 1
                    init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a + a/2, (j + 1/3)*a, (k + 1/3)*a + a/2])))
                    N += 1
                    init_pos = np.vstack((init_pos, np.array([(i + 1/3)*a, (j + 1/3)*a + a/2, (k + 1/3)*a + a/2])))
                    N += 1

    return init_pos, N, box_len


def init_velocity(N, D, T):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    N : int
        The number of particles in the system.
    D : int
        The number of spatial dimensinos of the system
    T : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel : np.ndarray
        Array of particle velocities
    """
    # center and standard deviation of gaussian distribution
    gauss_center = 0
    st_dev = T**0.5

    # generate velocities and substract mean for each dimension
    vel = np.random.normal(gauss_center, st_dev, (N,D))
    mean_vel = np.mean(vel, 0)
    vel -= mean_vel

    return vel


def rescale_velocity(N, T, t, num_tsteps, ekin, vel):
    """
    Function that rescales the velocities to correspond to a certain temperature
    Parameters
    ----------
    N : int
        Number of particles
    T : int
        The (dimensionless) temperature
    t : int
        The current timestep
    num_tsteps : int 
        Number of timesteps
    ekin : float
        Current kinetic energy of the system
    vel : (N x D) np.ndarray
        The current velocities of the particles

    Returns
    -------
    vel : (N x D) np.ndarray
        The rescaled velocities of the particles
    """
    goal_ekin = (N-1) * 3/2 * T
    ekin_error = goal_ekin / 20

    # conditions:
    # time = between 0 and 25% every 5%
    # difference ekin-goal > error
    if t > 0 and t <= (3*num_tsteps / 4) and t % (num_tsteps / 20) == 0 and np.absolute(ekin - goal_ekin) >= ekin_error:

        # scaling factor
        speed = np.sqrt(np.sum(np.power(vel, 2), 1))
        speed_square_sum = np.sum(np.power(speed, 2))
        scalefactor = np.sqrt(((N-1)*3*T/speed_square_sum))

        vel = vel*scalefactor
        print('rescaled\n')

        return vel, True

    return vel, False