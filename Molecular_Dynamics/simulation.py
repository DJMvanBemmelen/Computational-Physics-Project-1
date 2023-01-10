"""

This module can be considered the "back end" of the time evolution of the
simulation. It contains all functions that handle calculations needed for
the time evolution of the system and computing and saving data about the
observables of the system. The central function containing the main time
loop is simulate(). When the specified number of timesteps have been
simulated, simulate() saves the following data as .npy files:

- kinetic energy
- potential energy
- virial (for pressure calculation)
- relative distance histogram values (for pair correlation calculation)
- mean-squared displacement

A separate directory is created for each simulation. All of these directories
can be found in the folder 'Simulation Data'.

"""

import numpy as np
import math as m
from scipy.optimize import curve_fit
import os
import time

import initial_conditions as ic


def simulate(N, number_unitcells, D, T, pos, vel, box_len, dens, num_tsteps, tstep, bin_size, method, name):
    """
    Handling the time evolution of the system. This function uses either the
    Euler method or the Verlet method to update the position and velocity array's,
    of N interacting particles and collects data for computing the observables.

    Parameters
    ----------
    N : int
        Number of particles
    number_unitcells : int
        Number of unitcells of the center faced lattice
    D : int
        Number of dimensions
    T : int
        Dimensionless temperature of the system
    pos : (N x D) np.ndarray
        The initial positions of the atoms in Cartesian space
    vel : (N x D) np.ndarray
        The initial velocities of the atoms in Cartesian
    box_len : float
        The dimension of the simulation box
    dens : float
        The density of the system
    num_tsteps : int
        Number of timesteps
    bin_size : float
        Size of the histogram bins used in the pair correlation calculation
    tstep : float
        Duration of a single simulation step
    method : string
        Defines the numerical method (Euler or Verlet) 
    name : string
        Sets name of the directory for the outputted data

    Returns
    -------
    None
    """
    # time simulation start
    print('\nSimulation Progress:')
    start = time.time()

    # Initializing empty arrays for plotting quantities
    t_array = np.zeros(num_tsteps, dtype=float)
    ekin_t = np.zeros(num_tsteps, dtype=float)
    epot_t = np.zeros(num_tsteps, dtype=float)
    virial_t = np.zeros(num_tsteps, dtype=float)
    msd_t = np.zeros(num_tsteps, dtype=float)

    # Setting input for pair correlation
    n_bins = m.ceil(m.sqrt(D/4)*box_len/bin_size)
    rel_dist_hist_t = np.zeros((num_tsteps, n_bins), dtype=float)
    ref_particle = 0

    # Finding initial relative distances, positions and forces
    init_pos = pos
    rel_pos, rel_dist = atomic_distances(N, D, pos, box_len)
    f_res = force_LJ(N, D, rel_pos, rel_dist)

    n_eq_tsteps = 0

    # Simulating time-evolution
    for t in range(num_tsteps):

        # Checking simulation progress
        if t % (num_tsteps / 4) == 0:
            print(int(t/num_tsteps * 100), '%')

        # Calculating energy, virial, pair correlation and msd
        ekin_t[t] = kinetic_energy(vel)
        epot_t[t] = potential_energy(N, rel_dist)
        virial_t[t] = virial(N, rel_dist)
        rel_dist_hist_t[t] = pair_correlation(ref_particle, N, D, rel_dist, bin_size, box_len, n_bins)
        msd_t[t] = mean_squared_displacement(N, init_pos, pos)

        # Rescaling velocity
        vel, rescaled = ic.rescale_velocity(N, T, t, num_tsteps, ekin_t[t], vel)
        if rescaled:
            n_eq_tsteps = t

        # updating time, positions, relative distances, forces & velocities
        t_array[t] = t*tstep

        # ----- EULER METHOD ----- #
        if method == "Euler":
            pos, rel_dist = Euler(N, D, box_len, tstep, pos, vel)

        # ----- VERLET METHOD -----#
        if method == "Verlet":
            pos, f_res, rel_dist = Verlet(N, D, box_len, tstep, pos, vel, f_res)

    # Determining runtime
    finish = time.time()
    runtime = round(finish - start, 2)
    print('100 % \n\nFinished  -  Runtime:', runtime, 'seconds\n')

    # Making data directory
    name_folder = "Simulation Data"
    name_dir = "../" + name_folder + "/" + name

    while os.path.exists(name_dir):
        name_dir = name_dir + " copy"

    os.makedirs(name_dir)

    # Saving data to data directory
    np.save(os.path.join(name_dir, 'kinetic_energy'), ekin_t)
    np.save(os.path.join(name_dir, 'potential_energy'), epot_t)
    np.save(os.path.join(name_dir, 'virial'), virial_t)
    np.save(os.path.join(name_dir, 'distance_histogram'), rel_dist_hist_t)
    np.save(os.path.join(name_dir, 'mean_squared_displacement'), msd_t)

    if D == 2:
        np.save(os.path.join(name_dir, 'final_pos'), pos)

    # Saving Input data to file
    with open(os.path.join(name_dir, 'input.txt'), "w") as file:
        file.write("===== The input data of the simulation ===== \n")
        file.write("Number of spatial dimensions:\n")
        file.write(f"{D}\n")
        file.write("Box length: [sigma (3.405 Angstrom)]\n")
        file.write(f"{box_len}\n")
        file.write("Density: [sigma^-3 (3.405 Angstrom)]\n")
        file.write(f"{dens}\n")
        file.write("Number of unit cells: \n")
        file.write(f"{number_unitcells}\n")
        file.write("Number of particles: \n")
        file.write(f"{N}\n")
        file.write("Number of time steps: \n")
        file.write(f"{num_tsteps}\n")
        file.write("Time step size: [sigma * sqrt(mass / epsilon)] [\n")
        file.write(f"{tstep}\n")
        file.write("Temperature [epsilon/konst_Boltzman]: \n")
        file.write(f"{T}\n")
        file.write("Number of equilibriated time steps: \n")
        file.write(f"{n_eq_tsteps}\n")
        file.write("Bin size for the pair correlation: \n")
        file.write(f"{bin_size}\n")
        file.write("Run time: \n")
        file.write(f"{runtime}\n")

    return


def atomic_distances(N, D, pos, box_len):
    """
    Calculate the relative position vectors and relative distance between all atoms.

    Parameters
    ----------
    N : int
        Number of particles
    D : int
        Number of dimensions
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_len : float
        The dimension of the simulation box

    Returns
    -------
    rel_pos : (N x N x D) np.ndarray
        Relative positions of particles
    rel_dist : (N x N) np.ndarray
        The distance between particles (The diagonal elements are 0 because these
        represent the distance between a particle and itself)
    """
    rel_pos = np.zeros((N, N, D), dtype=float)
    rel_dist = np.zeros((N, N), dtype=float)

    # relative position vectors incorporating periodic BC's
    for i in range(N):
        rel_pos[i] = (pos - pos[i] + box_len/2) % box_len - box_len/2

    rel_dist = np.sqrt(np.sum(np.power(rel_pos, 2), 2))

    return (rel_pos, rel_dist)


def force_LJ(N, D, rel_pos, rel_dist):
    """
    Calculates the force acting on each particle as a result of the potential 
    created by all other particles using the Lenard Jones Potential. 

    Parameters
    ----------
    N : int
        Number of particles
    D : int
        Number of dimensions
    rel_pos : (N x N x D) np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : (N x N) np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    f_res : (N x D) np.ndarray
        The net force acting on particle i due to all other particles

    """
    f_res = np.zeros((N, D), dtype=float)

    # FULL FORMULA: -dU/dr * (1/r) * r_vec
    # f_abs_i = -dU/dr * (1/r)
    term1 = np.divide(12, np.power(rel_dist, 14), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)
    term2 = np.divide(6, np.power(rel_dist, 8), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)
    f_abs_i = -4 * (term1 - term2)

    f_i = np.einsum('ij,ijk->ijk', f_abs_i, rel_pos)

    f_res = np.sum(f_i, axis=1)

    return f_res


def Euler(N, D, box_len, tstep, pos, vel):
    """
    Uses the Euler method and the function atomic_distances() to update the
    position and velocity array's and the force.

    Parameters
    ----------
    N : int
        Number of particles
    D : int
        Number of dimensions
    box_len : float
        The dimension of the simulation box
    tstep : float
        Duration of a single simulation step
    pos : (N x D) np.ndarray
        The positions of the particles in cartesian space as calculated by the previous iteration
    vel : (N x D) np.ndarray
        The velocity of the particles in cartesian space as calculated by the previous iteration

    Returns
    -------
    pos : (N x D) np.ndarray
        The new positions of the particles in cartesian space
    rel_dist : (N x N) np.ndarray
        Relative particle distances to be used for energy calculations
    """
    # Update the position with periodic BC's
    pos = (pos + vel*tstep) % box_len

    # calculating new positions, distances and force
    rel_pos, rel_dist = atomic_distances(N, D, pos, box_len)
    f_res = force_LJ(N, D, rel_pos, rel_dist)

    # Updating the velocity
    vel += f_res * tstep

    return pos, rel_dist


def Verlet(N, D, box_len, tstep, pos, vel, f_res):
    """
    Uses the Euler method and the function atomic_distances() to update the
    position and velocity array's and the force.

    Parameters
    ----------
    N : int
        Number of particles
    D : int
        Number of dimensions
    box_len : float
        The dimension of the simulation box
    tstep : float
        Duration of a single simulation step
    pos : (N x D) np.ndarray
        The positions of the particles in cartesian space as calculated by the previous iteration
    vel : (N x D) np.ndarray
        The velocity of the particles in cartesian space as calculated by the previous iteration
    f_res : (N x D) np.ndarray
        The net force acting on particle i due to all other particles as calculated by the previous iteration

    Returns
    -------
    pos : (N x D) np.ndarray
        The new positions of the particles in cartesian space
    rel_dist : (N x N) np.ndarray
        Relative particle distances to be used for energy calculations
    f_res : (N x D) np.ndarray
        The net force acting on particle i due to all other particles to be used in the next iteration
    """
    # Position x(t+h) using: x(t), v(t) and F(x(t))
    pos = (pos + vel*tstep + f_res*(tstep**2 / 2)) % box_len

    # Force F(x(t+h)) using x(t+h)
    rel_pos, rel_dist = atomic_distances(N, D, pos, box_len)
    f_res_new = force_LJ(N, D, rel_pos, rel_dist)

    # Velocity v(t+h) using: v(t), F(x(t)) and F(x(t+h))
    vel += (f_res + f_res_new) * (tstep / 2)

    return pos, f_res_new, rel_dist


def kinetic_energy(vel):
    """
    Calculates the total kinetic energy of the system

    Parameters
    ----------
    vel : (N x D) np.ndarray
        Velocity of particle

    Returns
    -------
    ekin_tot : float
        The total kinetic energy of the system.
    """
    # calculate the magnitude of the velocities
    speed = np.sqrt(np.sum(np.power(vel, 2), 1))

    # calculate total kinetic energy using E = 0.5v^2
    ekin_tot = 0.5 * np.sum(np.power(speed, 2))

    return ekin_tot


def potential_energy(N, rel_dist):
    """
    Calculates the total potential energy of the atomic system using the Lennard Jones potential

    Parameters
    ----------
    N : int
        Number of particles
    rel_dist : (N x N) np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    epot_tot : float
        The total potential energy of the system.
    """
    A = np.divide(1, np.power(rel_dist, 12), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)
    B = np.divide(1, np.power(rel_dist, 6), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)

    epot_tot = 0.5 * np.sum(4 * (A - B))

    return epot_tot


def virial(N, rel_dist):
    """
    Calculates the non averaged virial term <0.5 Sum(ij) r_ij dU_ij/dr_ij>.

    Parameters
    ----------
    N : int
        Number of particles
    rel_dist : (N x N) np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    virial:  float
        Value for the sum of al r_ij dU_ij/dr_ij

    """
    A = np.divide(-12, np.power(rel_dist, 12), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)
    B = np.divide(-6, np.power(rel_dist, 6), out=np.zeros((N, N), dtype=float), where=rel_dist!=0)

    # One of the 0.5* is due to the double count of r_ij and r_ji
    vir = 0.5 * 0.5 * np.sum(4 * (A - B))

    return vir


def pressure(N, T, T_std, virial_t, n_eq_tsteps):
    """
    Calculates the Pressure of the system

    Parameters
    ----------
    N : int
        Number of particles
    T : float
        Average temperature of the system
    T_std : float
        STD of the Average temperature
    virial_t : np.ndarray
        The virial of the system at avery time iteration
    num_eq_tsteps : int
        Number of timesteps needed to let the system equilibriate

    Returns
    -------
    P : float
        Pressure of the system in units density/epsilon
    P_std : float
        Error margin on the pressure
    """
    # Calculating the avarage over simulation time from equilibrium time onwards
    viral_average = np.mean(virial_t[n_eq_tsteps :])

    # Determining the standard diviation (STD) for <virial>
    tmax = 2000
    virial_average_std = error_estimation(virial_t[n_eq_tsteps :], tmax)

    # Calculating the pressure
    P = T - 1/(3 * N) * viral_average
    P_std = T_std + 1/(3 * N) * virial_average_std

    return P, P_std


def pair_correlation(i, N, D, rel_dist, bin_size, box_len, n_bins):
    """
    Makes a histogram of the relative distance of the particles

    Parameters
    ----------
    i : int
        The refrence particle
    N : int
        Number of particles
    D : int
        Number of dimensions
    rel_dist : (N x N) np.ndarray
        Relative particle distances as obtained from atomic_distances
    bin_size : float
        distance between two histogram edges dr
    box_len : float
        The dimension of the simulation box
    n_bins : 

    Returns
    -------
    rel_dist_hist :  np.ndarray
        Number of particle pairs that have a distance between [r, r + dr]
    edges : np.ndarray
        List of distance tickes of the histogram edges [r_i, r_i + dr]
    """
    rel_dist_hist, edges = np.histogram(rel_dist[i], bins=n_bins, range=(0, m.sqrt(D/4)*box_len))

    # Remove the "distance" between the particle and itself
    rel_dist_hist[0] -= 1

    return rel_dist_hist


def pair_correlation_function(N, D, rel_dist_hist_t, n_eq_tsteps, bin_size, box_len):
    """
    Determines the pair correlation function g(x) that gives the number of pairs with a certain distance r

    Parameters
    ----------
    N : int
        Number of particles
    D : int
        Number of dimensions
    rel_dist_hist_t : (t_steps x n_bins) np.ndarray
        Histogram of the number of particle pairs that have a distance between [r, r + dr]
    num_eq_tsteps : int
        Number of timesteps needed to let the system equilibriate
    bin_size : float
        distance between two histogram edges dr
    box_len : float
        The dimension of the simulation box

    Returns
    -------
    g_r :  np.ndarray
        Pair correlation function that gives number of particle pairs that have a distance between [r, r + dr]
    edges : np.ndarray
        List of distance tickes of the histogram edges [r_i, r_i + dr]
    """
    # Calculating the avarage pair distance histogram over simulation time from equilibrium time onwards
    rel_dist_hist_average = np.mean(rel_dist_hist_t[n_eq_tsteps :], axis=0)

    # Determine error for each bin
    n_bins = m.ceil(m.sqrt(D/4)*box_len/bin_size)

    # Determine edges of bins for r
    edges = np.arange(0, m.sqrt(D/4)*box_len + m.sqrt(D/4)*box_len/n_bins, m.sqrt(D/4)*box_len/n_bins)
    r = edges + 0.5 * bin_size
    r = np.delete(r, -1)
    if len(r) == len(rel_dist_hist_average) + 1:
        r = np.delete(r, -1)
    V = box_len**D
    g_r = (2 * V)/(N * (N - 1)) * 1 / (4 * m.pi * bin_size) * np.divide(rel_dist_hist_average, np.power(r, 2))

    
    return r, g_r


def specific_heat(N, ekin_t, n_eq_tsteps):
    """
    Calculates specific heat using the Lebowitz theorem.

    Parameters
    ----------
    N : int
        Number of particles
    ekin_t : (1 x N) np.ndarray
        Total kinetic energy over simulation time
    num_eq_tsteps : int
        Number of timesteps needed to let the system equilibriate

    Returns
    -------
    C_v : float
        Specific heat of the system
    c_v : float
        Specific heat per particle
    """
    # Calculating the average over simulation time from equilibrium time onwards
    K_average = np.mean(ekin_t[n_eq_tsteps :])
    K2_average = np.mean(np.power(ekin_t[n_eq_tsteps :], 2))

    # Determining the standard diviation (STD)
    tmax = 2000

    # STD of <K^2>
    K2_average_std = error_estimation(np.power(ekin_t[n_eq_tsteps :], 2), tmax)

    # STD of <K>^2
    K_average_std = error_estimation(ekin_t[n_eq_tsteps :], tmax)
    K_average_2 = K_average**2
    K_average_2_std = K_average_2 * m.sqrt(2 * K_average_std/K_average)

    # Error of A =  <K^2> / <K>^2
    A = K2_average/K_average_2
    A_std = A * m.sqrt(K2_average_std/K2_average + K_average_2_std/K_average_2)

    # Calculating total specific heat
    C_v = 1/(1 + 1/N - A)
    C_v_std = C_v * A_std / A

    # Calculating specific heat per particle
    c_v = C_v/N
    c_v_std = C_v_std/N

    return C_v, C_v_std, c_v, c_v_std


def mean_squared_displacement(N, init_pos, pos):
    """
    Calculates the mean squared displacement.

    Parameters
    ----------
    N : int
        The number of particles in the system
    init_pos : N x D np.ndarray
        The position of each particle at the first timestep
    pos : N x D np.ndarray
        The position of each particle at the current timestep

    Returns
    -------
    msd : float
        The mean-squared displacement at the current timestep
    """
    dist_vec = pos - init_pos
    dist_squared = np.sum(np.power(dist_vec, 2), 1)

    msd = (1/N) * np.sum(dist_squared)

    return msd


def error_estimation(data, tmax):
    """
    Determines the standard deviation of the expectation value of the mean of a correlated data set

    Parameters
    ----------
    data : (1 x N) np.ndarray
        The data set containing N datapoints
    tmax : int
        The maximum value of t for determining the autocorrelation time.
        Relevant because for higher values t the autocorrelation function
        starts to fluctuate.

    Returns
    -------
    std : float
        The standard deviation of the expectation value of the mean of the data set
    """
    N = data.size
    t_array = np.arange(0, tmax, 1)
    autocor_func = np.zeros(tmax, dtype=float)

    # Determining the autocorrelation function
    for t in t_array:
        if t == 0:
            An, Ant = data, data
        else:
            An, Ant = data[0:-t], data[t:]

        numerator = (N-t)*np.sum(An*Ant) - np.sum(An)*np.sum(Ant)
        denominator = np.sqrt((N-t)*np.sum(np.power(An, 2)) - np.power(np.sum(An), 2)) * np.sqrt((N-t)*np.sum(np.power(Ant, 2)) - np.power(np.sum(Ant), 2))
        autocor_func[t] = numerator / denominator

    # Exponential function for fitting
    def exponential_function(t_array, autocor_time):
        return np.exp(-t_array/autocor_time)

    # Determining the autocorrelation time and standard deviation
    autocor_time, uncertainty = curve_fit(exponential_function, t_array, autocor_func)

    mean = np.sum(data)/N
    mean_squares = np.sum(data**2)/N
    std_dev = float(np.sqrt(2/N * autocor_time * (mean_squares - mean**2)))

    return std_dev


def temp_determination(N, ekin_t, n_eq_tsteps):
    """
    Determination of the average temperature after the rescaling stops

    Parameters
    ----------
    N : int 
        Number of particles
    ekin_t : (1 x N) np.ndarray
        Total kinetic energy over simulation time
    num_eq_tsteps : int
        Number of timesteps needed to let the system equilibriate

    Returns
    -------
    T_average : float
        Average temperatur of the system 
    T_average_std : float
        Standard deviation of the Temperature average
    """

    ekin_average = np.mean(ekin_t[n_eq_tsteps :])
    tmax = 2000
    ekin_average_std = error_estimation(ekin_t[n_eq_tsteps :], tmax)

    T_average = 2/3 * 1/(N - 1) * ekin_average
    T_average_std = 2/3 * 1/(N - 1) * ekin_average_std

    return T_average, T_average_std