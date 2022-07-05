#!/bin/bash
# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=percolate
#SBATCH --time=0:20:0
#SBATCH --exclusive
#SBATCH --cpus-per-task=1

# Replace [budget] below with your unique budget code
#SBATCH --account=
#SBATCH --partition=standard
#SBATCH --qos=standard

# Load the default HPE MPI environment
module load intel-compilers-19
module load mpt

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

srun -n 1 ./percolate
srun -n 2 ./percolate
srun -n 4 ./percolate
srun -n 8 ./percolate
srun -n 16 ./percolate
srun -n 32 ./percolate



