#!/bin/bash
#SBATCH -c 4                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 4-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=200G                        # Memory total in MB (for all cores)
#SBATCH -o roaryouterr/%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e roaryouterr/%j.err                 # File to which STDERR will be written, including job ID 
roary -f ./roary/$@ -e -n -p 8 -v gffs_sep/$@/*.gff