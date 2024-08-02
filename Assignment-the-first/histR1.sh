#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --output=logs/R1_%j.out 
#SBATCH --error=logs/R1_%j.err


conda activate bgmp_py312

./histogram.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz  -l 101 -o R1_hist.png


exit

