import pandas as pd
import numpy as np
import os

## NOTE - GEMMA EXECUTABLE PATH FOR FARHAT LAB!
gemma_executable = "/n/data1/hms/dbmi/farhat/SuperDCA/GEMMA/gemma-0.98.1-linux-static"
level = 'genedrop'
specieslists = [['E.coli', 'S.enterica'], ['M.tb', 'S.aureus']]
# iterate over antibiotics, species in each antibiotic
for idx, anti in enumerate(['cip', 'rif']): 
    for species in specieslists[idx]:
        namestem = f"{anti}_{level}_{species}"
        loci = namestem + ".loci"
        pheno = namestem + ".phenotypes"
        outfile = f"{anti}_{level}_{species}_gemma_out"
        grmname = f"{anti}_{level}_{species}.grm.txt"
        # run the actual gwas
        string = f"{gemma_executable} -g {loci} -p {pheno} -k {grmname} -maf 0 -lmm -notsnp -o {outfile}"
        os.system(string)
