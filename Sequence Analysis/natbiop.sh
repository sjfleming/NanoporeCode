#!/bin/bash

#SBATCH -J natbiop
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p general
#SBATCH --mem-per-cpu 1000
#SBATCH -t 1-00:00
#SBATCH -o Out/parout.txt
#SBATCH -e Out/parerr.txt

matlab -nodisplay -nosplash -r "natbio_params();exit"
