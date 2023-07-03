#!/bin/bash
#SBATCH --mail-user=akumar@nsstc.uah.edu  
#SBATCH -J wrf_plots
#SBATCH -p shared
#SBATCH --ntasks 8
#SBATCH -t 20-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=END,FAIL
#SBATCH -o slurm-%J.out # STDOUT
#SBATCH -e slurm-%J.err # STDERR

########## Add all commands below here

export PATH="/rhome/akumar/anaconda3/envs/pywrf/bin:$PATH"
#python plot_rainfall_houston.py
#python plot_vertical_cross_section_actual_vert_wind.py
#python plot_vertical_cross_section_quiver.py
#python plot_vertical_cross_section_ws.py
python plot_vertical_cross_section_actual_vert_wind_contourf.py
echo '**********Finished!**********'

