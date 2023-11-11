#!/bin/bash
#SBATCH --job-name=ribodesigner_parallel
#SBATCH --partition=scavenge
#SBATCH --time=04:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=jcm10@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=0-49

module load GCCcore/8.3.0
module load Python/3.7.2

## this lacks the installation of the pairwise analysis binary or packages
## so that should be taken into account

srun python3 main.py $1 $2 $3 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT