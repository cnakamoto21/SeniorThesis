import pandas as pd
import numpy as np
import os

gemma_executable = "/n/data1/hms/dbmi/farhat/SuperDCA/GEMMA/gemma-0.98.1-linux-static"
level = 'inter'
species = 'all'
for anti in ['cip', 'rif']: # 
    for analysistype in ['dropsite', 'classic', 'genedrop']: 
        namestem = f"{anti}_{level}_{species}"
        loci = namestem + ".loci"
        pheno = namestem + ".phenotypes"
        outfile = f"{anti}_{analysistype}_{species}_gemma_out"
        grmname = f"{anti}_{analysistype}_{species}.grm.txt"
        covname = namestem + ".covariates"
        # run the actual gwas
        string = f"{gemma_executable} -g {loci} -p {pheno} -k {grmname} -c {covname} -maf 0 -lmm -notsnp -o {outfile}"
        os.system(string)
