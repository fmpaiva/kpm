#!/bin/bash

#SBATCH --time=1:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core

module load gcc/gcc-14.2.0
g++ -O3 -mavx -I eigen-5.0.1/ src/*.cpp -o out
./out