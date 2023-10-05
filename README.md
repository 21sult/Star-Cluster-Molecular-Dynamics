# Star-Cluster-Molecular-Dynamics
Research on star cluster dissolution using molecular dynamics. The computational part is done in C++, and the data analysis is done in Python with Matplotlib. It is highly suggested that one read the research paper .PDF's if they wish to understand the context behind the code, as the following descriptions only give brief outlines.

## Description
In this work, I simulate the dynamics of a star cluster as a function of the distribution of masses and velocities.

A cluster is a gravitationally bound star system. In the absence of external forces, its structure is defined by the balance between the gravitational interaction that tends to bring the stars together and the distribution of random velocities that tend to disperse them. For sufficiently long times, the distribution of the speed of the stars must follow a Maxwell-Boltzmann distribution, as explained by energy equipartition theorem. Thus, lower mass stars should, on average, have superior velocities. When these speeds exceed the escape velocity, the respective stars can cease to be part of the cluster. As they escape, the mass of the cluster decreases, which leads to a decrease in escape velocity, allowing stars with slightly larger masses to also be able to leave the cluster, gradually leading to total dissolution of the system.

## Simulation of the interaction between two stars
Consider two stars, with masses m1 and m2, respectively. If we consider that stars are material particles, the interaction between them can be described as a simple gravitational attraction.
Here I Implement Verlet's method to simulate the interaction between two stars of unit solar mass. Initially, it was considered that the two are at a distance of (2 pc), have an initial speed of equal intensity (2 × 10**-2 km/s), opposite directions and direction perpendicular to the axis that passes through the center of the two. I obtain graphs of kinetic energy, potential energy and mechanical energy as a function of time.

## Simulation of the dynamics of the cluster
Attached you will find three files corresponding to initial conditions for a stellar cluster of N = {10,100,500} stars. Each file contains seven columns (x, y, z; vx, vy, vz; mass). For each case, I calculated the total number of stars at a distance less than r_t from the center of the cluster as a function of time for a time interval of ten million years.

## Generation of the initial configuration of the star cluster for the simulation
When formed, stars in a stellar system follow a mass distribution known as Initial Mass Function (IMF), which can be approximated by the Kroupa distribution:

![image](https://github.com/21sult/Computational-Physics/assets/145617965/c93b4780-a8bf-49fe-bbbd-a1770ecaf83b)

Where the value of alpha depends on the range of mass values of the star (in units of solar mass), with $\alpha = 0.3$ for $m < 0.08$, $\alpha = 1.3$ for $0.08 < m < 0.5$ and $\alpha = 2.3$ for $m > 0.5$. I started by generating a set of stars that follow this mass distribution. For simplicity, I assumed that the mass of one star is limited: $0.01 < m < 100$.

For the spatial distribution of stars, consider King's profile,

![image](https://github.com/21sult/Computational-Physics/assets/145617965/1d352dee-f192-44f3-bdd8-d6c004336687)

Where $n(r)$ is the number of stars at distance $r$ from the center of the cluster, $n_0$ is a scaling factor, $r_c$ is the radius of the core ($r_c = 1$ pc) and $r_t$ is the radius from which the density is zero ($r_t = 5$ pc).

The velocity distribution follows the Maxwell-Boltzmann distribution,

![image](https://github.com/21sult/Computational-Physics/assets/145617965/8c835b7e-1a05-45d2-a1d0-3aa8d420227d)

Where $v_p$ establishes the typical speed scale. I generated the speed of each star, assuming that $v_p = 1$ km/s. Please note that the direction of the velocity vector is uniformly distributed (https://mathworld.wolfram.com/SpherePointPicking.html). Calculated the average over several initial configurations for 100 stars. See the .pdf for all the findings.
