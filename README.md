# Bimono_Becker_Doring
Simulation of the bi-monomeric Becker-Döring system

# Prerequisite
* Matlab and methods about curves intersections [InterX](https://fr.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections).

# Description
This code is designed to simulate a bi-monomeric Becker-Döring system. 
The system is a large Ordinary Differential Equations (ODE) system to simulate the evolution in time of the concentration of the clusters and the monomers.
The system can be approximated by an advection-diffusion Partial Differential Equation (PDE) for the concentration of the clusters considering the size as a continuous variable coupled to two ODE for the concentration of monomers.
The mathematical properties of the system are described in [DFMR](https://www.sciencedirect.com/science/article/pii/S0022519319303194?casa_token=DCigPl-yfxQAAAAA:TjqSebXOEMsljjXYqvs2kri0VnzMaSxv6rZVdf3Xs41WxL9HeLp_YuLRpdXLiZseGqxMFyttaJGk) and in [DFMV](https://hal.science/hal-04635659).

# Numerical integration
The algorithms are based on :
*  Finite difference method for the PDE-ODE system using the Thomas algorithm [tridiagonal](https://github.com/tamaskis/tridiagonal-MATLAB),
*  Runge-Kutta 8th order integration for the ODE system.

# Running the tests
Save all the files in a repository and select the repository when you set the path in Matlab. Either run the file "all_phases_discrete.m" for the ODE system or the file "phase1_continuous.m" for the PDE/ODE system (it can take some time to run).
The code returns a series of graphs and an animation of the evolution of the concentration of the monomers in the phase space, the evolution of the "energy" and the evolution of the size distribution of the clusters.
