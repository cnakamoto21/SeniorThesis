## CODE TO MOVE THE NEW FAA FILES INTO THE FOLDER
import os
import subprocess as sp
import pandas as pd

# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
gids = list(set(list(patric_table.genomeid)))
for gid in gids:
    # make the copy
    cmd = 'cp ./prokkas/{n}/{n}.gff ./gff_files'.format(n = gid)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    print(out_cmd.decode())
