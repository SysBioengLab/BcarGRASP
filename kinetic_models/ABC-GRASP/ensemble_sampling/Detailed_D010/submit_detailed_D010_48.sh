#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=BcarDetailedD010n48
# Archivo de salida
#SBATCH --output=output_detailed_D010_48.txt
# Partición (Cola de trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --nodelist=n2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mail-user=benjamin.elizondo@uc.cl
#SBATCH --mail-type=ALL

export LD_PRELOAD=$HOME/anaconda3/lib/libstdc++.so.6
module load MATLAB/R2021a
srun matlab -nodisplay < b_car_rejection_detailed_D010_run.m