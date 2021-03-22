## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd

def grm_from_aln(input_data, idlists, species = "M. tuberculosis"):
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
    None; binary variant numpy file is produced
    
    """    
    
    ## MAKE BINARY-ENCODED MATRICES OF GENNOMES
    # read in alignments
    genomes = np.load("g_file.npy")
    genomes = np.transpose(genomes)
    # sort isolates in alphabetical order by PATRIC genome id
    PATRIC = pd.read_csv('carter_total_ids_mapped_to_rollingdb.csv', dtype = str)
    PATRIC = PATRIC.sort_values(by=['genome.genome_id'])['genome.biosample_accession']
    # read in biosample numbers from lab database identifier file
    BS = pd.read_csv('isolate_indexes.csv', dtype = str)
    BS = list(BS.reset_index(drop=True)['isolate_ID'])
    # get ordered list of indices in the variant matrix
    indices = [BS.index(asc) for asc in PATRIC]
    # index the matrix so that the order is correct
    genomes = genomes[indices]
    # pull out the reference
    snp_annotation = pd.read_pickle(snp_annotation_file)
    refseq = np.array([nuc_to_num[nuc] for nuc in snp_annotation['ref']])
    # make a new matrix of major and minor allele binaries
    bi_genomes = np.zeros((genomes.shape[0], genomes.shape[1]))
    # iterate over the rows, initiate binary matrix
    for i in range(bi_genomes.shape[0]):
        bi_genomes[i,:] = (refseq != genomes[i,:]) * 1
        
    # save the binary alignment numpy file
    np.save(f"grms/{species}_bi.npy", bi_genomes)

# set inputs for function calling: dictionaries, data files
nuc_to_num = {"A":0, "C":1, "G":2, "T":3, "N":9}
num_to_nuc = {0:"A", 1:"C", 2:"G", 3:"T", 9:"N"}
input_data = "g_file.npy"
patric_table = pd.read_csv('final_combined_genomes.csv', dtype = str)
idlist = list(set(list(patric_table.genomeid)))
idlist.sort()
snp_annotation_file = "/n/data1/hms/dbmi/farhat/Roger/homoplasy_project/rolling_DB_scrape/Genotypes_Filtered_2/genotypes_SNP_annotation.pkl"
# call function
grm_from_aln(input_data, idlist)
