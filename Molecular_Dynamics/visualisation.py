"""

This module contains all functions that are concerned with visualising
the data collected from simulations.

In main() users are to specify the name of the simulation to be visualised.
This name specifies the location of the data so it is important to check
the spelling to exactly match the name of the subfolder in Simulation Data.
The function input_reader(name) is called to read the input parameters from
input.txt. Main() also calls data_reader(name, <data name>) to read data
from the files and return it as np.arrays used by the plotting functions.

Main() is to be altered in order to create the wanted graphs.

"""

import numpy as np
import matplotlib.pyplot as plt
import os

import simulation as sim

def main():
    # Name of the simulation you want to visualise
    name = 'First_full_simulation'

    # Create path to directory
    name_folder = "Simulation Data"
    name_dir = "../" + name_folder + "/" + name

    # Getting input data from input.txt
    input_values, t_array = input_reader(name_dir)
    D = input_values[0]
    box_len = input_values[1]
    cel_len = input_values[2]
    number_unitcells = input_values[3]
    N = input_values[4]
    num_tsteps = input_values[5]
    tstep = input_values[6]
    T = input_values[7]
    n_eq_tsteps = input_values[8]
    bin_size = input_values[9]


    # Make plots
    # plot_energy(t_array, name_dir)
    plot_pair_correlation(name_dir, N, D, n_eq_tsteps, bin_size, box_len)
    # plot_msd(t_array, name_dir)

    # Get pressure and Specific heat info
    # pressure_info(name_dir, N, n_eq_tsteps)
    # specific_heat_info(name_dir, N, n_eq_tsteps)

    return


def pressure_info(name_dir, N, n_eq_tsteps):
    """
    Determines the pressure of the system

    Parameters
    ----------
    name_dir : string
        Name of directory where data is stored

    Returns
    -------
    The pressure of the system 
    """
    # Read data from .npy files
    virial_t = data_reader(name_dir, 'virial')
    P, P_std = sim.pressure(N, virial_t, n_eq_tsteps)

    print("The pressure of the system is equal to: ", P, "with an error of: ", P_std)
    return


def specific_heat_info(name_dir, N, n_eq_tsteps):
    """
    Determines the specific heat of the system 

    Parameters
    ----------
    name_dir : string
        Name of directory where data is stored

    Returns
    -------
    The pressure of the system 
    """

    # Read data from .npy files
    ekin_t = data_reader(name_dir, 'kinetic_energy')
    C_v, C_v_std, c_v, c_v_std = sim.specific_heat(N, ekin_t, n_eq_tsteps)

    print("The total specific heat of the system is equal to: ", C_v, "with an error of: ", C_v_std)
    print("The specific heat per particle is equal to: ", c_v, "with an error of: ", c_v_std)

    return


def plot_energy(t, name_dir):
    """
    Sums up the potential and kinetic energy and plots all three in one figure

    Parameters
    ----------
    t : np.ndarray
        Time stamps of every iteration
    ekin : np.ndarray
        Kinetic energy at every iteration
    epot : np.ndarray
        Potential energy at every iteration
    name_dir : string
        Name of directory where data is stored

    Returns
    -------
    Plots of the energy of the system over time.
    """
    # Read data from .npy files
    ekin = data_reader(name_dir, 'kinetic_energy')
    epot = data_reader(name_dir, 'potential_energy')

    # etot_array = np.zeros(num_tsteps, dtype=float)
    etotal = ekin + epot

    plt.plot(t, etotal, label="Total Energy", color='g')
    plt.plot(t, epot, label="Potential Energy", color='r')
    plt.plot(t, ekin, label="Kinetic Energy", color='b')
    plt.xlabel(r'Time [$\sigma\sqrt{\frac{m}{\epsilon}}]$')
    plt.ylabel(r'Energy [$\epsilon$]')
    plt.title('Energy')
    plt.legend()

    plt.savefig(os.path.join(name_dir, 'Energy'))
    plt.clf()

    return


def plot_atomic_distance2D(N, box_len, init_pos, name_dir):
    """
    Plots the initial condition with the relative distance vector

    Parameters
    ----------
    N : int
        Number of particles
    box_len : float
        The dimension of the simulation box
    init_pos : (N x 2) np.ndarray
        The initial positions of the atoms in 2D Cartesian space

    Returns
    -------
    Plots of the initial condition.
    """

    shift = np.array([-box_len, -box_len])
    shift = np.vstack((shift, np.array([-box_len, 0])))
    shift = np.vstack((shift, np.array([-box_len, box_len])))
    shift = np.vstack((shift, np.array([0, box_len])))
    shift = np.vstack((shift, np.array([box_len, box_len])))
    shift = np.vstack((shift, np.array([box_len, 0])))
    shift = np.vstack((shift, np.array([box_len, -box_len])))
    shift = np.vstack((shift, np.array([0, -box_len])))

    pos_mirrord = np.add(init_pos, shift[0])
    for i in np.arange(1, 8, 1):
        pos_mirrord = np.append(pos_mirrord, np.add(init_pos, shift[i]), axis=0)

    plt.clf()
    x, y = init_pos.T
    plt.scatter(x, y, color="b")

    x_m, y_m = pos_mirrord.T
    plt.scatter(x_m, y_m, color="y")

    plt.xlim(-box_len, 2 * box_len)
    plt.ylim(-box_len, 2 * box_len)
    plt.vlines(x=[0, box_len], ymin=-box_len, ymax=2 * box_len, color="k", linestyles="--")
    plt.hlines(y=[0, box_len], xmin=-box_len, xmax=2 * box_len, colors="k", linestyles="--")

    plt.savefig(os.path.join(name_dir, "Initial_position"))
    plt.clf()

    return

    
def plot_speed_distribution(vel, name_dir):
    """
    Plots the initial speed distribution of the particles. Used for checking
    if the distribution follows Maxwell-Boltzmann
    """
    speed = np.sqrt(np.sum(np.power(vel, 2), 1))

    n, bins, patches = plt.hist(speed, 100, density=True)

    plt.xlabel("speed")
    plt.title('100.000 particles at T = 300')

    plt.savefig(os.path.join(name_dir, "Velocity_distribution"))
    plt.clf()

    return


def input_reader(name_dir):
    """
    Reads the input data from the input.txt file

    Parameters
    ----------
    name : string
        Name of the directory containing the data files for the simulation

    Returns
    -------
    input_values : 1 x 8 np.ndarray
        The input values corresponding to the simulation
    """
    with open(name_dir + "/input.txt", "r") as input_file:
        lines = input_file.readlines()
        D = int(lines[2])
        box_len = float(lines[4])
        dens = float(lines[6])
        number_unitcells = int(lines[8])
        N = int(lines[10])
        num_tsteps = int(lines[12])
        tstep = float(lines[14])
        T = float(lines[16])
        n_eq_tsteps = int(lines[18])
        bin_size = float(lines[20])

    input_values = [D, box_len, dens, number_unitcells, N, num_tsteps, tstep, T, n_eq_tsteps, bin_size]

    t_array = np.linspace(0, num_tsteps*tstep, num_tsteps)

    return input_values, t_array


def data_reader(name_dir, data):
    """
    Read the data from the .npy file and returns as array for plotting

    Parameters
    ----------
    name_dir : string
        Name of the simulation repository that contains the data.npy files
    data : string
        Name of the data to be returned, without .npy

    Returns
    -------
    data_array : np.ndarray
        The requested data as array
    """
    data_array = np.load(os.path.join(name_dir, data + '.npy'))

    return data_array


if __name__ == '__main__':
    main()