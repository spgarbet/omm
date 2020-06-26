#!/bin/bash

#SBATCH --mail-user=jonny.schild@vumc.org
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --output=status/job%a.out
#SBATCH --error=status/job%a.out
#SBATCH --array=1-2
#SBATCH --constrain=haswell

[[ -d status ]] || mkdir status
[[ -d output ]] || mkdir output

module purge
module load Intel/2017.4.196
module load IntelMPI/2017.3.196
module load R/3.4.3

Rscript --no-save sim-accre.R ${SLURM_ARRAY_TASK_ID}

