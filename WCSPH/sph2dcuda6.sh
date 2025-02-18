#!/bin/bash
#SBATCH --job-name=SPH2D
#SBATCH --partition=unlimited
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --mail-type=end
#SBATCH --mail-user=schakr18@umd.edu
module load cuda
module load gcc
module load openmpi
./SPH2DCUDA ./SPHinputFileInitialization.txt