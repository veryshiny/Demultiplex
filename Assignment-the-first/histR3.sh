#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --output=logs/R3_%j.out 
#SBATCH --error=logs/R3_%j.err

conda activate bgmp_py312

./histogram.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz  -l 8 -o R3_hist.png


exit

