## CODE TO DOWNLOAD ALL OF THE FNA FILES   
import os
import pandas as pd
import subprocess as sp

# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv')
gids = list(set(list(patric_table.genomeid)))
for genome_name in gids:
    # command line search
    cmd = 'wget -O fafiles/{n}.fa -qN ftp://ftp.patricbrc.org/genomes/{n}/{n}.fna'.format(n = str(genome_name))
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    if idx % 100 == 0:
        print("completed", idx)
