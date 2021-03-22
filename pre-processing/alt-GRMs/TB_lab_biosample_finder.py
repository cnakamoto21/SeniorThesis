## identify isolates that overlap between my dataset and lab dataset
## code from Anna Green, using data / code from other lab members

import pandas as pd
# Read table of all isolates in RollingDB
# NB replace /Users/agreen/Farhat_O2_mount with /n/data1/hms/dbmi/farhat to run on o2
all_isolates = pd.read_csv("/n/data1/hms/dbmi/farhat/anna/genome_logreg/input_data/isolate_annotation.csv", index_col=0)
# Read table of patric isolates from Carter
patric_isolate = pd.read_csv("combined_genomes.csv",dtype=str)
# filter to just tb
tblab = [i.split('.')[0] == "1773" or i.split('.')[0] == "1733" for i in list(patric_isolate.genomeid)]
patric_isolate = patric_isolate.loc[tblab]
### Read in all the PATRIC genomes and their biosample_accessions
patric_table = pd.read_csv("/n/data1/hms/dbmi/farhat/anna/patric_data/MTB_list.csv", sep="\t", dtype=str)
patric_table['id'] = patric_table['genome.genome_id']
# Carter's table with only isolates that we can map to a biosample
patric_m = patric_isolate.merge(patric_table, left_on='genomeid', right_on='id', how="left")
patric_m = patric_m.dropna(subset=["genome.biosample_accession"], axis=0)
patric_m["biosample"] = patric_m["genome.biosample_accession"]
# make an output table
ids_in_rollingdb = set(all_isolates.isolate_ID)
patric_m = patric_m.query("biosample in @ids_in_rollingdb")
patric_m.to_csv("carter_total_ids_mapped_to_rollingdb.csv")
