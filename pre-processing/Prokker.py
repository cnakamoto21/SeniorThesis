# Actual shell code
import os
import pandas as pd
import subprocess as sp

# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
gentoid = {"562": 'Excherichia', "1773": 'Mycobacterium', "1733": 'Mycobacterium', "28901": 'Salmonella', "1280": 'Staphylococcus'}
gids = list(set(list(patric_table.genomeid)))
for gn in gids:
    # get the name of the genome
    genus = gentoid[gn.split('.')[0]]
    cmd = f'sbatch run_prokka.sh {gn} {genus}'
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    if idx % 100 == 0:
        print("completed", idx)
