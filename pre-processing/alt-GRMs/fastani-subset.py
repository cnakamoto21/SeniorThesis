# Actual shell code
import os
import pandas as pd
import subprocess as sp
import sys

if __name__ == "__main__":
    # read in the table downloaded from patric
    patric_table = pd.read_csv('adj_combined_genomes.csv', dtype = str)
    idx = int(sys.argv[1])
    row = patric_table.iloc[idx]
    gn1 = row.genomeid
    # subset by species
    species_ind = [en.split('.')[0] == row.genomeid.split('.')[0] for en in patric_table['genomeid']]
    patric_table = patric_table.loc[species_ind]
    for idx2, row2 in patric_table.iterrows():
        # skip if comparison already done
        if idx2 <= idx:
            continue
        gn2 = row2.genomeid
        if os.path.exists(f"fastanis/{gn1}_{gn2}.out"):
            continue
        # call fastani
        cmd = f'FastANI-master/fastANI -q fna_files/{gn2}.gff -r fna_files/{gn1}.gff -o fastanis/{gn1}_{gn2}.out'
        p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
        out_cmd, err = p.communicate()
    if idx % 100 == 0:
        print("completed", idx)    
