#!/bin/bash
#SBATCH --job-name=ribodesigner_parallel
#SBATCH --partition=scavenge
#SBATCH --time=04:00:00
#SBATCH --export=ALL
#SBATCH --mail-user=jcm10@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=0-9

module load GCCcore/12.3.0
module load Python/3.11.3
pip install --user -r requirements.txt

## this lacks the installation of the pairwise analysis binary or packages
## so that should be taken into account

srun python3 main.py $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT