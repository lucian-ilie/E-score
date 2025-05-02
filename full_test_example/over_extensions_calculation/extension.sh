#!/bin/bash

#SBATCH --job-name=EXTENSION        # Job name
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (no gpus needed)
#SBATCH --mem=16G                   # Memory per nodes
#SBATCH --time=01:30:00             # Time limit hrs:min:sec

srun ./extensions.py "cd01068" "MODEL"

