#!/bin/bash 


# >>> select --> Per node
# >>> ncpus --> Core per node
# >>>  mpiprocs --> Core per node
#PBS -l select=25:ncpus=1 -l place=free

# >>> Job name 
##PBS -N wlcN

# >>> Queue name 
# >>> Defalut queue --> workq
#PBS -q workq

#PBS -r n 

#PBS -j oe 

# Load environments.
module purge
module load intel/2021.4.0

NPROC=`wc -l < $PBS_NODEFILE`

cd $PBS_O_WORKDIR

mpicc turrouse.c -o turrouse -lm
mpicc msxrouse.c -o msxrouse -lm
mpicc longturrouse.c -o longturrouse -lm
mpicc longmsxrouse.c -o longmsxrouse -lm
mpicc shortturrouse.c -o shortturrouse -lm
mpicc shortmsxrouse.c -o shortmsxrouse -lm


for f in 10 50 
do
for gammar in 10 
do
for v in 0 10 50 100
do
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n $NPROC   ./turrouse $N $v $l $dt $gammar $f > output.txt
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n 1        ./msxrouse $N $v $l $dt $gammar $f > output.txt
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n $NPROC   ./shortturrouse $N $v $l  $gammar $f > output.txt
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n 1        ./shortmsxrouse $N $v $l  $gammar $f > output.txt
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n $NPROC   ./longturrouse $N $v $l 0.1 $gammar $f > output.txt
mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n 1        ./longmsxrouse $N $v $l 0.1 $gammar $f > output.txt
done
done
done

#mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n $NPROC ./shortturrouse $N $v $l  $gammar $f > output.txt
#mpirun -bootstrap ssh -machinefile $PBS_NODEFILE -n 1      ./shortmsxrouse $N $v $l  $gammar $f > output.txt

