#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --output=logs/R4_%j.out 
#SBATCH --error=logs/R4_%j.err


conda activate bgmp_py312

./histogram.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz  -l 101 -o R4_hist.png


exit

