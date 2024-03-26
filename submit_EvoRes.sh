#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=50G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --job-name=EK_EvoRes

module load singularity/3.5.2

singularity run evores_0.2.1.sif "5 - DistributionComparisons.R"
