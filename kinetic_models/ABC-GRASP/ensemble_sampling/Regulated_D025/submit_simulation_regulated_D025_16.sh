#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=SimulationRegulatedD025n16
# Archivo de salida
#SBATCH --output=simulation_regulated_D025_16.txt
# Partición (Cola de trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --nodelist=n2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=benjamin.elizondo@uc.cl
#SBATCH --mail-type=ALL

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/anaconda3/lib
echo $LD_LIBRARY_PATH
module load MATLAB/R2021a
srun matlab -nodisplay < b_car_simulation_regulated_D025_run.m