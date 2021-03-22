## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import csv

# remake file with all genomes and only subset for tb based on what Anna had
genomepd = pd.read_csv('combined_genomes.csv', dtype = str)
tbpd = pd.read_csv('carter_total_ids_mapped_to_rollingdb.csv', dtype = str)
keeps = [i for i in genomepd.genomeid if i.split('.')[0] != "1773" and i.split('.')[0] != "1733"]
keeps += [i for i in genomepd.genomeid if i in tbpd.genomeid.values]
genomepd = genomepd.loc[genomepd['genomeid'].isin(keeps)]

## clean dataset based on alignment drops
with open("S. aureus_out_of_aln.txt", "r") as infile:
    drop_from_staph = infile.readlines()
drop_from_staph = [line.rstrip('\n') for line in drop_from_staph]
with open("S. enterica_out_of_aln.txt", "r") as infile:
    drop_from_sal = infile.readlines()
drop_from_sal = [line.rstrip('\n') for line in drop_from_sal]
droplist = drop_from_staph + drop_from_sal
genomepd = genomepd[~genomepd["genomeid"].isin(droplist)]

# get rid of the psky ass singly e coli
ec_ids = [i for i in genomepd.genomeid if i.split('.')[0] == "562"]
ec = genomepd.loc[(genomepd['genomeid'].isin(ec_ids)) & (genomepd['antibiotic'] == 'rifampicin')]
genomepd = genomepd.drop([ec.index[0]])

# save final new genome table for whole genome analyses
genomepd.to_csv('final_combined_genomes.csv', index=False)
