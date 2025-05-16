# TUR_AOUP

This code and data validate the TUR relation in mutli-particle active Ornstein-Uhlenbeck system.
See paper: https://arxiv.org/abs/2410.22126

# Code explain
- `submit_mpi.sh`

Bash script to initiate the MPI job to compile and run next 6 c codes. You can change the arguments of each code in this file.

ex) `qsub -N turr -v N=50,v=10,l=2,dt=0.01,gammar=10.0,f=10.0 submit_mpi.sh`

- `turrouse.c`, `longturrouse.c`, `shortturrouse.c`

c code to simulate the AOUP attached to the rouse chain with length of $2N+1$ particles. 
It gives the files `dump~.xyz` of trajectory of every particles every certain time steps,
and the files `taetaxdot~.xyz` of average value of $\langle x\circ \dot \eta \rangle$.
For single AOUP system, set `N=0` in this code. 
Prefix of `long~` or `short~` is for the long time simulation or short time simulation with different interation time.

- `msxrouse.c`, `longmsxrouse.c`, `shortmsxrouse.c`

Calcuate MSD from the trajectory ensembles generated from the previous code. 
It gives `msx~.dat`.


