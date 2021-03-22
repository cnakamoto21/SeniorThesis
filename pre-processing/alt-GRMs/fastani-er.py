# Actual shell code
import os
import pandas as pd
import subprocess as sp
# read in the table downloaded from patric
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
for idx, row in patric_table.iterrows():
    if idx < 4700:
        continue
    cmd = f'sbatch fastani-er.sh {idx}'
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    out_cmd, err = p.communicate()
    if idx % 100 == 0:
        print("completed", idx)    
    if idx > 10000:
        break
