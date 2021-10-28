#!/bin/sh
#SBATCH --job-name=test_assemblers
#SBATCH --output=slurm_stdout_stderr.txt
#SBATCH --ntasks=80
#SBATCH --mem=132000
#SBATCH --time=168:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends

python3 test_assemblers.py