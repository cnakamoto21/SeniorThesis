## IMPORT RELEVANT PACKAGES / FUNCTIONS
import numpy as np
import pandas as pd
import os
import subprocess as sp

genomepd = pd.read_csv('combined_genomes.csv', dtype = str)
ec_ids = [i for i in genomepd.genomeid if i.split('.')[0] == "562"]
genomepd = genomepd.loc[(genomepd['genomeid'].isin(ec_ids))]
GRM = np.empty((genomepd.shape[0],genomepd.shape[0]))
for idx, row in genomepd.iterrows():
    for idx2, row2 in genomepd.iterrows():
        if idx2 < idx:
            continue
        if idx2 == idx:
            GRM[idx,idx] = 1
            continue
        # get the name of the genome
        gn1 = row.genomeid
        gn2 = row2.genomeid
        if not os.path.exists(f"fastanis/{gn1}_{gn2}.out"):
            cmd = f'FastANI-master/fastANI -q fna_files/{gn2}.gff -r fna_files/{gn1}.gff -o fastanis/{gn1}_{gn2}.out'
            p = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
            out_cmd, err = p.communicate()
        with open(f"fastanis/{gn1}_{gn2}.out", "r") as infile:
            text = infile.read()
        GRM[idx,idx2] = float(text.split()[2]) / 100
        GRM[idx2,idx] = GRM[idx,idx2]

np.save(f"grms/E. coli_raw_grm.npy", GRM)
