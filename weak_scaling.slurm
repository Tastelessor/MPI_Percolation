#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=percolate
#SBATCH --time=0:20:0
#SBATCH --exclusive
#SBATCH --nodes=1

# Replace [budget] below with your unique budget code
#SBATCH --account=

# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard

# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Load the default HPE MPI environment
module load intel-compilers-19
module load mpt

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

srun -n 1 ./percolate 8759 0.4064 30
echo "-------------------------------"
srun -n 1 ./percolate 8759 0.4064 60
srun -n 2 ./percolate 8759 0.4064 60
echo "-------------------------------"
srun -n 1 ./percolate 8759 0.4064 120
srun -n 4 ./percolate 8759 0.4064 120
echo "-------------------------------"
srun -n 1 ./percolate 8759 0.4064 240
srun -n 8 ./percolate 8759 0.4064 240
echo "-------------------------------"
srun -n 1 ./percolate 8759 0.4064 480
srun -n 16 ./percolate 8759 0.4064 480
echo "-------------------------------"
srun -n 1 ./percolate 8759 0.4064 960
srun -n 32 ./percolate 8759 0.4064 960


