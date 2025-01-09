#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=Enzymes2DetailedD025n16
# Archivo de salida
#SBATCH --output=enzymes2_detailed_D025_16.txt
# Partición (Cola de trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --nodelist=n2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=benjamin.elizondo@uc.cl
#SBATCH --mail-type=ALL

export LD_PRELOAD=$HOME/anaconda3/lib/libstdc++.so.6
module load MATLAB/R2021a
srun matlab -nodisplay < b_car_simulation_enzymes_second_detailed_D025_run.m