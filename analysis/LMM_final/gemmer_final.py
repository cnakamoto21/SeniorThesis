import pandas as pd
import numpy as np
import os

## NOTE - MANUALLY CHANGED ALL OF THE FILE NAMES TO HAVE NO SPACES!!
gemma_executable = "/n/data1/hms/dbmi/farhat/SuperDCA/GEMMA/gemma-0.98.1-linux-static"
for idx, anti in enumerate(['ciprofloxacin', 'rifampicin']):
    for species in [['all', 'E.coli', 'S.enterica'], ['all', 'M.tb', 'S.aureus']][idx]:
        namestem = f"{species}_{anti}"
        loci = namestem + ".loci"
        pheno = namestem + ".phenotypes"
        outfile = namestem + "_gemma_out"
        grmname = namestem + ".grm.txt"
        # run the actual gwas
        string = ""
        if species == 'all':
            covname = namestem + ".covariates"
            string = f"{gemma_executable} -g {loci} -p {pheno} -k {grmname} -c {covname} -maf 0 -lmm -notsnp -o {outfile}"
        else:
            string = f"{gemma_executable} -g {loci} -p {pheno} -k {grmname} -maf 0 -lmm -notsnp -o {outfile}"
        if species == 'S.aureus':
            os.system(string)
            break
        if species == 'all':
            covname = namestem + ".covariates"
            string = f"{gemma_executable} -g fp_{loci} -p fp_{pheno} -k fp_{grmname} -c fp_{covname} -maf 0 -lmm -notsnp -o fp_{outfile}"
        else:
            string = f"{gemma_executable} -g fp_{loci} -p fp_{pheno} -k fp_{grmname} -maf 0 -lmm -notsnp -o fp_{outfile}"
        os.system(string)
