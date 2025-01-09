#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MetabolitesDetailedD010n16
# Archivo de salida
#SBATCH --output=metabolites_detailed_D010_16.txt
# Partición (Cola de trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --nodelist=2
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=benjamin.elizondo@uc.cl
#SBATCH --mail-type=ALL

export LD_PRELOAD=$HOME/anaconda3/lib/libstdc++.so.6
module load MATLAB/R2021a
srun matlab -nodisplay < b_car_simulation_metabolites_detailed_D010_run.m