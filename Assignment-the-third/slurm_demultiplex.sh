#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 24
#SBATCH --mem=100G
#SBATCH --output=logs/demultiplex_%j.out 
#SBATCH --error=logs/demultiplex_%j.err


conda activate bgmp_py312

/usr/bin/time  -v ./script_demultiplex.py -i /projects/bgmp/shared/2017_sequencing/indexes.txt -f1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -f2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -f3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -f4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz 

exit