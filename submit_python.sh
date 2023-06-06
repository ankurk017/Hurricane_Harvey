#!/bin/bash
#SBATCH --mail-user=akumar@nsstc.uah.edu  
#SBATCH -J CRIDA
#SBATCH -p shared
#SBATCH --ntasks 8
#SBATCH -t 20-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=END,FAIL
#SBATCH -o slurm-%HRRR02.out # STDOUT
#SBATCH -e slurm-%HRRR02.err # STDERR

########## Add all commands below here

export PATH="/rhome/akumar/anaconda3/envs/pywrf/bin:$PATH"

python F1_S1_extract_rainfall.py
#python india_locations_finer_for_districts_april03.py
echo '**********Finished!**********'

