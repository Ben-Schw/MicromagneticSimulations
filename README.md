This repository contains projects about numerical solutions in statistical physics.
The main focus lies on problems in micromagnetics.

We have small projects surrounding Monte Carlo simulations of a 2D or 3D (anti)ferromagnet. It investigates the temperature dependence of the order parameter, which is also the output of a CSV file. (language: C)

Also, this repostiory contains a cluster analysis with the Hoshen-Kopelman algorithm. In a 2D lattice the sites are occupied at a certain probability. We compute the largest cluster depending on this probability. (language: C)

Furthermore, there is a folder containing smaller issues like computing Pi through Monte Carlo simulation.

Lastly, we have the LLG solver (Landau-Lifshitz-Gilbert equation). This simulates the dynamics of individual spins in a lattice. Global temperature, magnetic field, and uniaxial anisotropy is implemented. The grid assumes isotropic exchange between the spins, which is the same for all spins. The parameters of this simulation can be altered via "parameters.h".

Future add-ons might contain: 
- Implementation of sublattices of the grid
- Locality of magnetic and temperature field
- ...
