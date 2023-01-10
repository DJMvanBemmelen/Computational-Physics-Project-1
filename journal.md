# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Monday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
(due 14 February 2022, 23:59)

### Plan of Attack
1. To store each particles position and velocity at each time step we will create arrays: `pos`, `vel`. Both of size N x D where N is the number of particles and D is the number of spatial dimensions.

2. To calculate the force on the i'th particle, we first have to calculate the distances to all other particles. We'll write a function: `atomic_distances` that performs this calculation using the `pos` array as input. This function will take into account the periodic boundary conditions using function `periodic_bc_2D`. We will then wirte a function, `force_LJ`, that uses the output of `atomic_distances` and the gradient of the Lennard-Jones potential formula to generate the resulting force on the i'th particle by summing over the forces from all other particles. 

3. To calculate the new position and velocity of the i'th particle, the funtion 'Euler' method will implement the numerical method called the Euler method. This function will also take into account the periodic boundary conditions.

4. To calculate the total energy, we will write a function, `kinetic`, that calculates the total kinetic energy using `vel` as input. We'll also write a function, `potential_LJ`, that calls `atomic_distances` and uses the output to calculate the total potential energy in the system using the Lennard-Jones potential formula. Finally we'll find the total energy of the system by summing both energy's. 


### Progress 
1. Made as discussed in Plan of Attack. This can be seen at the beginning of [week1.py](https://hub.compphys.quantumtinkerer.tudelft.nl/user/dvanbemmelen/lab/tree/Project-1_mwortelboer/Weekly%20Progress/week1.py). It was however unclear to us if a function that generates the initial conditions in the form of a fcc latice position and gaussian velocity distibution needed to be made this week.

2. Made as discussed in Plan of Attack. In the figure the way the `atomic_distances` are calculated using the `periodic_bc_2D` is displayed. The way the force is calculated can be also be found in [week1.py](https://hub.compphys.quantumtinkerer.tudelft.nl/user/dvanbemmelen/lab/tree/Project-1_mwortelboer/Weekly%20Progress/week1.py). However, With the current dimensions the force becomes very small. The interactions between the particles are therefore not seen very clear. Also the way the `force_LJ` and `potential_LJ` are calculated can possibly be improved by removing the for loop. This is something we would like to discuss in class. 

![Atomic Distance](Figures/Atomic_distance.png) 

3. Made as discussed in Plan of Attack. See [week1.py](https://hub.compphys.quantumtinkerer.tudelft.nl/user/dvanbemmelen/lab/tree/Project-1_mwortelboer/Weekly%20Progress/week1.py)

4. Made as discussed in Plan of Attack. See point 2 for discussion on `potential_LJ`. See [week1.py](https://hub.compphys.quantumtinkerer.tudelft.nl/user/dvanbemmelen/lab/tree/Project-1_mwortelboer/Weekly%20Progress/week1.py)

5. (EXTRA) Finally all files are managed in a way that seemed appropiate in the folder [Interacting Particles Lennard Jones @?@](https://hub.compphys.quantumtinkerer.tudelft.nl/user/dvanbemmelen/lab/tree/Project-1_mwortelboer/Interacting%20Particles%20Lennard%20Jones%20%40%3F%40). However, we are not sure if this is correct. We would like to discuss this in class. The way constants can be made `global` is also something we need to find out still. Lastly, we were also woundering if it would be better if we would make the symulations via a class instead of a function. This is the last thing we would like to discuss in class. 

## Week 2
(due 21 February 2022, 23:59)

### Plan of Attack
As can be seen in the progres of week1, there are a lot of specific week2 functionalitys already implemented. For example the minimal image convention, plotting posibility's for the energy aswel as the plotting functionality for the atomic distances. The rest of the Plan of Attack for week 2 wil therefor concern optimalisation of these functions and the implementation of the third demention. 

1. Fix the way the files are is structured. Fix the `global` variables. 

2. Change the functions: `force_LJ()`, `potential_LJ()` and `kinetic()` into their nondimensional counterparts. Also change the `global` variables into the nondimensional ones.

3. Maybe rewrite implementation of minimal image convention in `atomic_distances()` using the modulo operator.

4. Investigate implementation of the modulo operator in the time development of position to automatically implement periodic boundary conditions, making `interacting_particles.periodic_bc()` obsolete.

5. Do the proposed simulation and plot the calculated energy.

6. Remove the for loop from the function `force_LJ()`


### Progress 

1. The files are structured in a way that is more clear by dividing the functions into three files instead of 5. This can be seen in the folder [Molecular_Dynamics](Project-1_mwortelboer/Molecular_Dynamics). The functions: `atomic_distances()`, `force_LJ()`, `periodic_bc_2D()`, `periodic_bc_3D()`, `potential_LJ()`, `kinetic()` and `Euler()` are now all put into this file [interacting_particles.py](../interacting_particles.py) whereas the initial condition functions will from now on be clustered into the file [initial_condition.py](../initial_condition.py) and the different simulations will be clustered into the file [simulation.py](../simulation.py). Also the folder weekly progress is deleted and from now on all coding progress will be made in the [Molecular_Dynamics](Project-1_mwortelboer/Molecular_Dynamics) folder. 

2. We removed the variables `mass`, `sigma` and `epsilon` from the functions `force_LJ()`, `potential_LJ()` by setting them to 1. In the function `kinetic()` we did however keep an epsilon??? Furthermore we changed the `box_len` to be equal to the number of particles `N` and changed the timestep to be equal to `tstep = 1E-3`. 

3. By implementing the modulo operator in the way proposed in the lecture notes of week 2, we managed to significantly shorten the function atomic_distances() in the module interacting_particles.py. The new way to determine the relative positions of particles inherently accounts for the possibility of the particle being closer to an image across the boundary. Therefore the several lines of codes including nested loops we first used to perform this check are no longer needed. By further examining this function we got rid of another loop used to calculate the relative distance between particles, increasing the functions efficiency while reducing the lines of code needed. 

4. The implementation of the modulo operator in the function Euler() was very simple and very effective. By replacing `pos += vel*tstep` by `pos = (pos + vel*tstep) % box_len`, we've ensured the resulting position will always be between 0 and `box_len`. This makes the whole function periodic_bc() obsolete and is therefore a massive improvement of our code.

5. With al the changes to the code we managed to simulate the interacting Argon particles. Sadly however, the energy in the system is not conserved as can be seen in the first plot with `tstep = 1E-3` and the epsilon in the function `kinetic`. When we tried to change the the time step to `tstep = 1E-4` and we removed the epsilon from the function `kinetic` not much changed for the energy conversion. The weird thing is the plot of the potential, where it seems to go from zero to inf every time step. We understand that currently our code does not work properly. We believe that the problem either lays within the calculation of the force, or the implementation of the nondimensionalisation or the boundary conditions. We would like to have a look at this with a TA. 
    
![Energy_Plot_chaos](Figures/Energy_Plot_chaos.png) 

![Energy_Plot_epsilon_1](Figures/Energy_Plot_epsilon_1.png) 

6. By using the numpy function `np.einsum('ij,ijk->ijk',2D,3D)` for the operation where the prefactor $f_abs_i = - \frac_{dU}{dr} \frac{1}{r}$ is multiplied with the 3D array `rel_pos`. 


## Week 3
(due 28 February 2022, 23:59)

### Plan of Attack
There are some milestones this week that we've already implemented: Our code is already compatible to take on more than 2 particles and we have already decided on a structure for our repository that we feel is orderly and easy to read. We would like to have feedback on the structure if it is not perceived this way. Then there are still issues with the milestones of last week that will be incorporated in our plan of attack.

1. Identify and fix bug in code causing the unrealistic energy profiles. 
    - Together
2. Implement the velocity-Verlet algorithm.
    - Daniel
3. Investigate the energy profiles of the system and compare Verlet to Euler
    - Mees

### Progress

1. We have identified the bug in our system to be in `force_LJ()`. The multiplication between the prefactor of the force and the relative position vectors was incorrect, yielding an opposite sign of what was to be expected. This in turn led to an attractive force for particles very close together and a repulsive force for particles very far away, resulting in enormous accelerations and exploding energies. Luckily there is an easy fix for this. In our method the force prefactor from particle A on particle B is multiplied by the relative position vector from B to A, instead of the vector from A to B. S the size is correct but the direction is reversed. By simply adding a minus sign to the final step of the force calculation we've fixed our bug. In the image below the energies are plotted for 2 particles in a box of size (6*sigma)<sup>3</sup>. The particles are positioned such that they are just inside the repulsive range of the force, i.e. just left of the minimum of the LJ-potential curves. They are given a very small initial velocity toward eachother and the position and velocity are updated with the Euler Method. You can nicely see the initial slowing down of the particles to zero kinetic energy due to the repulsive force. Then the rapid repulsion and again slowing down, this time due to the the attractive force. We also see the periodic behaviour that we would expect for these initial conditions. We believe this verifies that our implementation of the force and the energies are correct.

![Energy_Euler_2P](Figures/Energy_Euler_2P.png)

2. Implemented the velocity-Verlet algorithm by calculating the force on t = 0 as can be seen in the file [simulation.py](../simulation.py) and then updating every new force within the function Verlet which can be found in the file [interacting_particles.py](../interacting_particles.py). Here the periodic boundary is implemented the same way as discussed in point 4 of week 2. 

3. By implementing the velocity-Verlet algorithm we can clearly see improvement in the energy conservation of our system. The images below show the energy curves of a system of 10 particles in a box of size (7*sigma)<sup>3</sup> give random initial positions and velocities. At first glace the difference between the Euler Method and the Verlet Method isn't very apparent. But looking at the bottom two pictures it becomes clear that the total energy variation of the Verlet Method is a factor 100 smaller. 

![Energy_Euler_10P](Figures/Energy_Euler_10P.png) ![Energy_Verlet_10P](Figures/Energy_Verlet_10P.png)

![Energy_Total_Euler_10P](Figures/Energy_Total_Euler_10P.png) ![Energy_Total_Verlet_10P](Figures/Energy_Total_Verlet_10P.png)

## Week 4
(due 7 March 2022, 23:59)

### Plan of Attack

1. Implement the initialization of positions onto an face centered cubic lattice in the function: `fcc_lattice()`. 

2. Implement the initialization of velocities obeying a Maxwell-Boltzmann distribution in the function: `init_velocity()`. 

3. Implement the rescaling of temperature at a certain timestep inside the function:`simulation()` by making a function `rescale_velocity()` in the file [initial_conditions.py](../initial_conditions.py). Furthermore we will show how the desired temperature is attained after a certain amount of rescaling and equilibrating by making energy plots.

4. Study at least one observable, and compare its behaviour to literature. We will decide on which one later in the week. 

### Progress

1. Implemented the fcc lattice for 2D and 3D where a 2D example can be seen in the next figure. There are a few things worth noting. Firstly the number of unit cells `number_unitcells` is given as an argument instead of the number of particles. Furthermore, the lattice constant `a` is derived from the box length and the number of unit cells. Both these decisions insure that the final result is an homogeneus distribution of the particles throughout the main box and its mirrored sides. Finally, to ensure that there are no particles on the boundary of the box all the particles are shifted by [1/3, 1/3, 1/3]. 

![Initial Position fcc](Figures/Initial_position_fcc.png)

2. We implemented the initialization of velocities making use of numpy's built in function `numpy.random.normal()` which generates an array of prescribed size where all numbers are drawn from a normalized gaussian distribution with inputted center and standard deviation. For the center we took 0, because particles are allowed to have negative velocities for the distinct spatial directions. For the standard deviation we took $sd = \sqrt{T}$. We non-dimensionalized T using the factor from the gaussian distribution, so T is now measured in units of $[\frac{\epsilon}{k_{B}}]$. From the lecture notes of week one we know $\frac{\epsilon}{k_{B}} = 119.8 [K]$ so room temperature corresponds to approximately $T = 2.5$. We substract the center-of-mass speed from all velocities by first calculating the average for each spactial direction with the `numpy.mean()` function and then substracting it from the velocity of each particle. In the image below you can see a histogram of the speeds ou function produced for 100.000 particles at T = 300. It appears to follow the Maxwell-Boltzmann distribution nicely

![velocity distribution](Figures/Velocity_distribution.png)

3. We implemented the rescaling of the velocities to match the inputted temperature in the function `rescale_velocity()`. Using our dimensionless T the formula for the goal kinetic energy from temperature becomes: $$v^2 = (N-1)3T$$. We call this function in `simulate()` in the file `simulation.py`. For the first quarter of the simulation it checks every 5% of the total nr of timesteps if the kinetic energy of the system is within 1 percent of the goal kinetic energy. This means our system can have a maximum of 5 rescalings.

![velocity rescale](Figures/Energy_rescale.png)

4. The observables Pressure and Pair correlation where implemented in the code. (i) For the Pressure the nondimensionalisation led to a pressure in terms of [$\epsilon \rho$] which nondimensionalised the full equation to $$P = 1 - \frac{1}{3N}\left< \sum_{ij}r_{ij}\frac{\partial U}{\partial r_{ij}} \right>.$$ This equation was implemented in the function: `virial_theorem()`. The term that is averages over simulation time $\sum_{ij}r_{ij}\frac{\partial U}{\partial r_{ij}}$ is calculated in the function: `virial()` by using the `rel_dist` for every timestep. (ii) For the Pair correlation the relative distances of all $N-1$ particles w.r.t particle 1 are stored in the array `rel_dist_hist_t`. Then the function: `pair_correlation_function()` averages the distance of every k'th particle w.r.t particle 1 over all equilibrium timesteps. It then counts the number of particles for which the distance lies between the pre-fixed $[r, r + \Delta r]$. This number then represents the $\left<n(r) \right>$. For a initial simulation with (N = 256, lattice distance a = 1.4 $[\sigma]$ and T = 3 $[\frac{\epsilon}{k_{b}}]$) we got values of P = 4.376 $[\rho \epsilon]$ and the a pair correlation value between [3,4] of $g(r) = 0.0094$

## Week 5
(due 14 March 2022, 23:59)

### Plan of Attack 
1. Implement calculation of errors in the function: `error_estimation()` and test your implementation on data with a known correlation time. You can get the code to generate random data with a specified autocorrelation time
- Mees
2. Compute observables including errors. Specifically we need to implement the observables: Diffusion and Specific heat 
- Daniel
3. Make sure your code is structured logically. (IDEA: Implement the printing of info to text files and implement the running of visualisations into another function ????) 

4. Make a plan for simulations to go into the report: How do you want to validate your simulation, and which observables/simulations do you want to run?

5. Make a brief statement about the efficiency of your code: can you run the simulations you would like to do? Would it make sense to improve performance?

### Progress
1. We chose to implement the autocorrelation function method for the error estimation. We verified our implementation in the file Error_estimation.ipynb on randomized data and the results were satisfactory. The only part of this that we're still uncertain about is the tmax value. In the lecture notes it is stated that computing the autocorrelation function for higher values of t, since the autocorrelation function has very large fluctuations for this range. We also saw this behaviour in the notebook. We are unsure what value for t is a good maximum and how we can determine this. Coming lecture we would like to get feedback on this.

2. We had implemented the pressure and the pair correlation function last week. For the pressure the implementation of the error was very straightforward and proposed no problems. The pair correlation error proposed more problems. We seem to get an error related to the tmax value discussed in the first point, but this error only starts to come up in the pair correlation function. We would also like to discuss this with a TA coming lecture. We implemented the specific heat and the error without problems. The mean-squared displacement was also pretty straight forward to compute and the error was implemented without big troubles. But for all of the observables the tmax value is still an uncertainty for us.

3. 
    - We finalized the structure of our repository and rewrote a lot of functions. All of our relevant code is now in the folder 'Molecular_Dynamics'. Here we have the files `simulation.py`, containing all functions needed for running the simulation such as the calculation of the force, the verlet method for time evolution, etc.; `visualisation.py` containing all functions for plotting; `observables.py`, containing all functions for calculating observables, like the energy, pair correlation, pressure, etc. and finally `initial_conditions.py`, containing the functions needed for initializing the system.
    - We also rewrote our code to improve the handling of data in our simulation. First we had to alter simulation.py to plot different observables and this took a lot of time rerunning simulations when we wanted a different observable, or we would have to compute all observables. To work around this we did the following: We have made a folder 'Simulation Data' in our directory. When a simulation is run a folder is made specifically for that simulation in 'Simulation Data'. All data for all observables and the input given to the simulation is saved in .npy files in this folder carrying the simulation name. Visualisation.py contains functions that can read out this data and transfer it back to numpy.arrays and uses this data to plot graphs. This extra step of saving and reading the data may seem like extra work, but it saves us a lot of time when trying to alter graphs, working on functions in visualisation.py. Furthermore it is good practice to save the raw data for later use. This alteration was made with the idea in mind that this simulation could be used after this project and to be a framework for future projects on a larger scale. 
    
4. The first observable we want to discuss in our report is the pressure in the system. We want to run multiple simulations keeping everything constant but the temperature and then plot the (T,P) curve. If our simulation works properly we should see pressure going up with increasing temperature. We also hope this might provide some inside in phase transistions, or at least regions of distinctly different phase, but we need to do more research in what we exactly expect to observe. The pair correlation function could also be a good measure for this. We also want to investigate the varying timesteps and see at what size the data starts to converge, i.e. at what point reducing the timestep no longer gives large improvement in the accuracy of the simulation. 

5. We made an effort to reduce the use of for-loops in favor of numpy array operations as much as possible. As python is simply a thin wrapper around C++, using numpy allows for so much more efficient computing that it results in a computation time only a factor 2 of using pure C++ code, but with way simpler and shorter code. In our opinion we've succeded quite well in this. Places where we didn't manage are for example in initializing the center faced chrystal lattice, but since this is done only once per simulation it does not act as a bottleneck for our computation time. As discussed in point 3. we also made an effort to improve the data handling of our simulation which also makes postprocessing data way more efficient. Our code runs a simulation of 10,000 time steps and a couple hundred particles in about a minute. We feel like further improving our computation efficiency would require great effort and would not amount to large improvements.