#!/bin/sh
#SBATCH --job-name=test_assemblers
#SBATCH --output=slurm_stdout_stderr.txt
#SBATCH --ntasks=80
#SBATCH --mem=132000

python3 test_assemblers.py