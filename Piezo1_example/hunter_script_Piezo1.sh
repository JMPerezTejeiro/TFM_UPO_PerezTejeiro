#! /usr/bin/env bash

#SBATCH --mem='60gb'
#SBATCH --constraint=cal
#SBATCH --time='10:00:00'
hostname

data_folder=~/hunter/Piezo1
results_folder=~/hunter/Piezo1/results_f06/

# Carga de librerías (desde cuenta de software)
source ~soft_bio_267/initializes/init_degenes_hunter

# 1. Análisis de expresión diferencial
degenes_Hunter.R -i $data_folder'/count_table.txt' -o $results_folder -f 0.6 -t target_file -m DE


# 2. Análisis funcional
functional_Hunter.R -i $results_folder -t E -c 6 -o $results_folder/Piezo1_functional_enrichment -m Mouse
