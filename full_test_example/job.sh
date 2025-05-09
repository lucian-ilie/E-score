#!/bin/bash

#SBATCH --job-name=cd01068         #Job name
#SBATCH --cpus-per-task=1          #Number of CPU cores per task
#SBATCH --gres=gpu:1               #Number of GPUs per task 
#SBATCH --mem=10G                  #Memory per nodes
#SBATCH --time=1:30:00            #Time limit hrs:min:sec


srun ./testing_code.py cd01068
