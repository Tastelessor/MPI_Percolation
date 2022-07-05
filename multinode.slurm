#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=percolate
#SBATCH --time=0:20:0
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --tasks-per-node=20
#SBATCH --cpus-per-task=1 
#SBATCH --ntasks=40

# Replace [budget] below with your unique budget code
# For EPCC Students: [budget]=dc134-[student-id]
# For External Students: [budget]=dc135-[student-id]
#SBATCH --account=dc134-s2168291

# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard

# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Load the default HPE MPI environment
module load intel-compilers-19
module load mpt

# Change to the submission directory
cd $SLURM_SUBMIT_DIR

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.

srun ./percolate 8759


