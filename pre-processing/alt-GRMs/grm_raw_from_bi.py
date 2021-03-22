## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import sys

def grm_from_bi(bi_mat, species = "S. aureus"):
    """
    This function makes a genetic relatedness matrix (GRM) out of the binary variant numpy matrix for one species. 
    
    """    
    
    ## MAKE BINARY-ENCODED MATRICES OF GENNOMES
    # read in binary matrix
    binary = np.load(bi_mat)
    grm = np.matmul(binary, np.transpose(binary)) / binary.shape[1]
        
    # save the grm
    np.save(f"grms/{species}_raw_grm.npy", grm)

species_list = ["S. aureus", "S. enterica", "M. tuberculosis"]

if __name__ == "__main__":
    i = int(sys.argv[1])
    grm_from_bi(f'grms/{species_list[i]}_bi.npy', species_list[i])
