## get isolates from lab data!!
## code from Anna Green, using data / code from other lab members

import numpy as np
import pandas as pd
from scipy import stats

# Input information from Roger
genotype_translator = {
0:"A", 1:"C", 2:"G", 3:"T", 9:"N"
}
isolate_annotation_file = "/n/data1/hms/dbmi/farhat/Roger/homoplasy_project/rolling_DB_scrape/Genotypes_Filtered_2/genotypes_isolate_annotation.pkl"
genotype_array_file = "/n/data1/hms/dbmi/farhat/Roger/homoplasy_project/rolling_DB_scrape/Genotypes_Filtered_2/genotypes_matrix.npy"
snp_annotation_file = "/n/data1/hms/dbmi/farhat/Roger/homoplasy_project/rolling_DB_scrape/Genotypes_Filtered_2/genotypes_SNP_annotation.pkl"

isolate_ids = pd.read_csv("../input_data/isolate_annotation.csv")

isolates_to_keep = pd.read_csv("../input_data/carter_ids_mapped_to_rollingdb.csv", sep=",")
isolates_to_keep = set(isolates_to_keep.biosample.values)

k = isolate_ids.query("isolate_ID in @isolates_to_keep")
k.to_csv("isolate_indexes.csv")

### Grabs only the columns corresponding to isolates in our dataset - needs a lot of memory so don't rerun unless needed
ga = np.load(genotype_array_file)
g = ga[:,k.index]
np.save("g_file.npy", g)
