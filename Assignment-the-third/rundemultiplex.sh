#!/bin/bash
#SBATCH --account=bgmp               ###
#SBATCH --partition=bgmp       ### Partition
#SBATCH --output=demultiplex%j.out
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --mail-type=END              ### Mail events (NONE, BEGIN, END, FA$
#SBATCH --mail-user=sobermil@uoregon.edu  ### Where to send mail
#SBATCH --cpus-per-task=1
#SBATCH --error=demultiplex%j.err

conda activate bgmp_py310
cd /projects/bgmp/sobermil/bioinfo/Bi622/Assignment-the-third
/usr/bin/time -v ./demultiplex.py \
 -read1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz\
  -read2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
  -index1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz\
   -index2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
   -output output -known_indexes known_indexes.txt