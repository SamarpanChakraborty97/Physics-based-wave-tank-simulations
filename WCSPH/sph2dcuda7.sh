#!/bin/bash
#SBATCH --job-name=SPH2D
#SBATCH --partition=gpup100
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=4
#SBATCH --time=6:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=schakr18@umd.edu
module load cuda
module load gcc
module load openmpi
./SPH2DCUDA ./SPHinputFileNew.txt