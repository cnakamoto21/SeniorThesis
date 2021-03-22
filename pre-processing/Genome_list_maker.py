import pandas as pd

## CODE TO COMBINE ALL THE TABLES INTO ONE TALBE OF GENOME IDS
# read in individual dataframes
frames = []
for i in [562, 77643, 28901, 1280]: 
    frames.append(pd.read_csv(f"amr_data_{i}.tsv", sep='\t', dtype = str))
# combine into single frame
big = pd.concat(frames, ignore_index=True)
# drop unnecessary columns
big = big[['genome.genome_id', 'genome_drug.antibiotic', 'genome_drug.resistant_phenotype']]
# rename columns for ease
big = big.rename(columns={"genome.genome_id": "genomeid", "genome_drug.antibiotic": "antibiotic", 
                          "genome_drug.resistant_phenotype": "resistance"})
# drop rows with antibiotics we don't care about (note: already screened for biological test of resistance)
big = big.loc[(big.antibiotic == 'ciproflox') | (big.antibiotic == 'ciprofloxacin') | 
              (big.antibiotic == 'rifampin') | (big.antibiotic == 'rifampicin')]
# standardize genome names
big.loc[big.antibiotic == 'ciproflox', 'antibiotic'] = 'ciprofloxacin'
big.loc[big.antibiotic == 'rifampin', 'antibiotic'] = 'rifampicin'
# save as csv for future use
big.to_csv('combined_genomes.csv', index = False)
