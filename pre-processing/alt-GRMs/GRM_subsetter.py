import pandas as pd
import numpy as np

## MAKE GRM FILES FROM GENERAL GRMS
level = 'wgs'
# spclists = [['all', 'E.coli', 'S.enterica'], ['all', 'M.tb', 'S.aureus']]
spclists = [['all', 'S.aureus']]
allspclists = [['S.aureus', 'M.tb', 'S.enterica', 'E.coli'], ['S.aureus', 'M.tb']]
genomes = pd.read_csv("final_combined_genomes.csv", dtype = str)
idlist = list(set(list(genomes.genomeid)))
nametosp = {"S.aureus": '1280', "E.coli": '562', "M.tb": '1773', "S.enterica": '28901'}

# function to subset the raw grm for each species based on ids in analysis of interest
def grm_subsetter(species, idlist, anti, gids):
    if species != "E.coli":
        idlist.sort()
    gidsub = [i for i in gids if i.split(".")[0] == nametosp[species] or (species == "M.tb") and i.split(".")[0] == '1733']
    tempid = [i for i in idlist if i.split(".")[0] == nametosp[species] or (species == "M.tb") and i.split(".")[0] == '1733']
    gididx = [tempid.index(gid) for gid in gidsub]
    grm = np.load(f'grms/{species}_raw_grm.npy')
    grm = grm[gididx][:,gididx]
    return grm
                
# iterate over antibiotics and species, calling 
for i, anti in enumerate(['ciprofloxacin', 'rifampicin']): 
    for species in spclists[i]:
        # read in gids in analysis
        with open(f"gemma_final/{species}_{anti}_gids.txt", "r") as infile: # ADD FP_ FOR FP
            gids = infile.readlines()
        gids = [line.rstrip('\n') for line in gids]
        # if all-species analysis, iterate over all species, calling subsetter and combining subsets
        if species == 'all':
            grm_subsets = []
            for sp in allspclists[i]:
                grm_subsets.append(grm_subsetter(sp, idlist, anti, gids))
            grm_n = sum([i.shape[0] for i in grm_subsets])
            grm = np.zeros((grm_n, grm_n))
            for idx, sub_grm in enumerate(grm_subsets):
                start = sum([i.shape[0] for i in grm_subsets[:idx]])
                grm[start:start+sub_grm.shape[0],start:start+sub_grm.shape[0]] = sub_grm
        # else just call the subsetter once
        else:
            grm = grm_subsetter(species, idlist, anti, gids)
        # save output in file format appropriate format
        np.savetxt(f"gemma_final/{species}_{anti}.grm.txt", grm, delimiter="\t") # ADD FP_ FOR FP
        print(len(gids), grm.shape)
