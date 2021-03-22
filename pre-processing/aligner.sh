#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-08:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40000                        # Memory total in MB (for all cores)
#SBATCH -o alnouterr/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e alnouterr/%j.err                 # File to which STDERR will be written, including job ID
./muscle -in $@.faa -out $@.afa -maxiters 2