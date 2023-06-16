#!/bin/bash
#SBATCH --mail-user=akumar@nsstc.uah.edu  
#SBATCH -J houston
#SBATCH -p standard
#SBATCH --ntasks 8
#SBATCH -t 20-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=END,FAIL
#SBATCH -o slurm-%HRRR02.out # STDOUT
#SBATCH -e slurm-%HRRR02.err # STDERR

########## Add all commands below here

export PATH="/rhome/akumar/anaconda3/envs/pywrf/bin:$PATH"

python plot_vertical_cross_section_quiver.py
#python india_locations_finer_for_districts_april03.py
echo '**********Finished!**********'

