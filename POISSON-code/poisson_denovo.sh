#!/bin/bash

#SBATCH -J denovo
#SBATCH -n 1
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu 2000
#SBATCH -t 0-24:00
#SBATCH -o Out/stdout_%a.txt
#SBATCH -e Out/stderr_%a.txt
#SBATCH --open-mode=append

matlab -nojvm -nodisplay -nosplash -r "poisson_seq(\${SLURM_ARRAY_TASK_ID});exit"
