#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4000                         # Memory total in MB (for all cores)
#SBATCH -o mouterr/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e mouterr/%j.err                 # File to which STDERR will be written, including job ID
/home/chn6/.conda/envs/Prokka_113/bin/prokka --kingdom Bacteria --outdir prokkas/$1 --genus $2 --prefix $1 fna_files/$1.fa