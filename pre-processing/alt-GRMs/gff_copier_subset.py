# Actual shell code
import os
import pandas as pd
import subprocess as sp

# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
gids = list(set(list(patric_table.genomeid)))
gentoid = {"562": 'Escherichia', "1773": 'Mycobacterium', "1733": 'Mycobacterium', "28901": 'Salmonella', "1280": 'Staphylococcus'}
for gn in gids:
    # get the name of the genome
    gn = row.genomeid
    genus = gentoid[gn.split('.')[0]]
    cmd = f'cp gff_files/{gn}.gff gffs_sep/{genus}'
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    if idx % 100 == 0:
        print("completed", idx)
