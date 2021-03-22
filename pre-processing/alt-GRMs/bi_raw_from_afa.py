## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd

def grm_from_aln(gen_aln, idlists, speciestosp, species = "S. aureus"):
    """
    This function makes a genetic relatedness matrix (GRM) out of a whole genome alignment for ONE species. 
    
    Parameters
    ----------
    gen_aln_names: list of strings
    The file location and names of the full genome alignment
    idlist: list of strings
    The gennome ids, in alphabetical order, that should be incorporated into the grm
    date: string
    The date of the run, in order to differentiate the product from others
    speciestosp: dict
    Converts from species name string (i.e. "M. tuberculosis") to taxon id in PATRIC (i.e. 1773)
    
    Returns
    -------
    None; binary variant matrix numpy file is produced
    
    """    
    
    ## MAKE BINARY-ENCODED MATRICES OF GENNOMES
    # read in alignments
    with open(gen_aln, "r") as infile:
        genomes = Alignment.from_file(infile, format="fasta")
    # pull out the reference 
    refseq = genomes[[idx for idx, name in enumerate(genomes.ids) if "ref" in name][0]]
    # subset and reorder alignment using premade idlist
    spidlist = [genomes.id_to_index[gid] for gid in idlist if gid.split(".")[0] == speciestosp[species] and gid in genomes.ids]
    genomes = genomes.select(sequences = spidlist)
    # make a new matrix of major and minor allele binaries
    bi_genomes = np.zeros((genomes.N, genomes.L))
    # iterate over the rows, initiate binary matrix
    for i in range(bi_genomes.shape[0]):
        bi_genomes[i,:] = (refseq != genomes.matrix[i,:]) * 1

    # get a list of the isolates out of the alignment to update dataset
    out_of_aln = [gid for gid in idlist if gid.split(".")[0] == speciestosp[species] and gid not in genomes.ids]
    with open(f"{species}_out_of_aln.txt", 'w') as f:
        for item in out_of_aln:
            f.write("%s\n" % item)
        
    # save the binary variant matrix numpy file
    np.save(f"grms/{species}_bi.npy", bi_genomes)

# set the inputs for the species with core genome alignments
speciestosp = {"S. aureus": "1280", "S. enterica": "28901"}
species_list = ["S. aureus", "S. enterica"]
gen_aln_list = ["roary/Staphylococcus/core_gene_alignment.aln", "roary/Salmonella/core_gene_alignment.aln"]
# get a list of the ids to try out
patric_table = pd.read_csv('combined_genomes.csv', dtype = str)
idlist = list(set(list(patric_table.genomeid)))
idlist.sort()
# iteratively call the function to make the binary variant matrix
for i in range(2):
    grm_from_aln(gen_aln_list[i], idlist, speciestosp, species_list[i])
