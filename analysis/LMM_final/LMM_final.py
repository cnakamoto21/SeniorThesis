## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import csv

def phen_cleaner(alignmentfile, rsstdata, antibiotic): 
    """
    This function generates cleaned phenotype data for a given antibiotic from a csv with the data
    
    Parameters
    ----------
    alignmentfile: string
    the file path and name of your alignment file (.afa)
    rsstdata: string
    the file path and name of a file with data on resistance (.csv)
    must have columns for genome ids (as in PATRIC), antibiotic, and resistance phenotype
    resistance phenotype should be string for each - "Susceptible", "Intermediate", "Resistant"
    antibiotic: string
    the name of the antibiotic of interest as it appears in the rsstdata file
    
    Returns
    -------
    list of integers
        indicates susceptible (0) or resistant (1) phenotypes in same order as numpy array rows. 
    matrix (1D) populated by bools
        each entry indicates whether or not the corresponding isolate has phenotype data for the antibiotic of interest
    """
    
    ## READ IN RELEVANT DATA
    # get an alignment object (with matrix and stuff) from alignment file
    with open(alignmentfile, "r") as infile:
        aln = Alignment.from_file(infile, format="fasta")
    # get all phenotype data
    genomepd = pd.read_csv(rsstdata, dtype = str)

    ## GET PHENOTYPE DATA, RESTRICT ACCORDINGLY
    # get the phenotypes in the same order as the matrix 
    phens = []
    keepers = np.empty((aln.matrix.shape[0]))
    gids = list(aln.ids)
    # iterate over genome ids in the alignment
    for i in range(aln.matrix.shape[0]): 
        # subset phenotype data based on whether or not antibiotic data at given id
        rsst = genomepd.loc[(genomepd.genomeid == gids[i]) & (genomepd.antibiotic == antibiotic)]
        # if no corresponding data, store index for removal later
        if rsst.resistance.values.shape[0] == 0:
            keepers[i] = False
            continue
        keepers[i] = True
        # otherwise get corresponding antibiotic resistance phenotype
        if rsst.resistance.values[0] == "Susceptible": 
            phens.append(0)
        else: 
            phens.append(1)
            
    ## RETURN THE PHENOTYPE DATA
    return phens, keepers

## function that generates clean binary variant matrices from alignment files
def aln_cleaner(alignmentfile, keepers = None, oldaln = None, 
                sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}):
    
    ## READ IN RELEVANT DATA
    # get an alignment object (with matrix and stuff) from alignment file
    with open(alignmentfile, "r") as infile:
        aln = Alignment.from_file(infile, format="fasta")    
    # if there is another alignment, reorder the current one to be the same and drop stuff not in that matrix
    alncln = aln
    if oldaln: 
        alncln = aln.select(sequences = [alncln.id_to_index[x] for x in oldaln.ids if x in alncln.ids])
    # if not, keepers applies! use it to drop bad indices
    else: 
        alncln = aln.select(sequences = keepers.astype(bool))
    
    ## CREATE BINARY MATRIX AND 90% SCREEN
    # create list of species for each isolate
    specieslabs = [i.split('.')[0] for i in list(alncln.ids)]
    # get reference sequences as an alignment object - hardcode out 1733 error
    seqs = [aln.id_to_index[f"{i}.ref"] for i in list(sptoid.keys()) if i != "1733"]
    refs = aln.select(sequences = seqs)
    tuple([alncln.identities_to(refs[i]) for i in range(4)]) # ran into error here
    np.vstack(tuple([alncln.identities_to(refs[i]) for i in range(4)]))
    # create a matrix to conduct the identity screen
    identities = np.vstack(tuple([alncln.identities_to(refs[i]) for i in range(4)])).T
    # make a new matrix of major and minor allele binaries
    muts = np.zeros((alncln.matrix.shape[0], alncln.matrix.shape[1]))
    muts.fill(np.nan)
    # iterate over the rows, initiate binary matrix and check identity level
    identityfails = []
    for i in range(muts.shape[0]):
        muts[i,:] = (refs.matrix[sptoid[specieslabs[i]]] != alncln.matrix[i,:]) * 1
        if identities[i,sptoid[specieslabs[i]]] < 0.1: # CHECK ON THIS - CHANGED IT FROM 0.9! SHOULD BE CORRECT NOW
            identityfails.append(alncln.ids[i])
    
    ## RETURN RELEVANT DATA
    return muts, alncln, refs, identityfails

## function that calls phen and aln cleaners to get cleaned combined binary variant matrices, phenotypes, and labels
def aln_combiner(alignmentfiles, rsstdata, antibiotic, 
                 sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}, 
                 sptoname = {"562": "E. coli", "1773": "M. tb", "1733": "M. tb", "28901": "S. enterica", "1280": "S. aureus"}):
    ## RUN OTHER CUSTOM FUNCTIONS TO GET THE ALIGNMENTS, PHENOTYPES
    # get phenotype data
    phens, keepers = phen_cleaner(alignmentfiles[0], rsstdata, antibiotic)
    mutlist, alnlist, reflist, fullidlist, allidfails = [], [], [], [], []
    # iterate over the genes, run the alignment cleaner for each
    for i in range(len(alignmentfiles)):
        if i == 0:
            muts, alncln, refs, idfails = aln_cleaner(alignmentfiles[i], keepers, sptoid = sptoid)
            fullidlist = alncln.ids
        else: 
            muts, alncln, refs, idfails = aln_cleaner(alignmentfiles[i], oldaln = alnlist[0], sptoid = sptoid)
            fullidlist = list(set(fullidlist).intersection(alncln.ids))
        # drop appropriate sites and add them to lists
        mutlist.append(muts)
        alnlist.append(alncln)
        reflist.append(refs)
        allidfails = allidfails + idfails

    ## RESTRICT ISOLATES IN MUTS, PHENOTYPES
    # get a list of ids that are beinng kept; reorder to be alphabetical
    idsinall = list(set(fullidlist).difference(set(allidfails)))
    idsinall.sort()
    # use list of kept ids to correct all
    for i in range(len(mutlist)):
        selection_index = [alnlist[i].id_to_index[x] for x in idsinall]
        mutlist[i] = mutlist[i][selection_index]
        alnlist[i] = alnlist[i].select(sequences = selection_index)
        # phenotype list ordered like first, so fix on that iteration
        if i == 0:
            phens = [j for i, j in enumerate(phens) if i in selection_index]

    # CREATE SPECIES LABELS
    splabs = [sptoname[i.split(".")[0]] for i in alnlist[0].ids]
    
    return mutlist, phens, splabs, reflist, idsinall

def fildrop(splabs, mutlist, phens, reflist, genenames, nametoid, anti, ids, species = None, 
            ref = "M. tb"):
    """
    This function makes the phenotype, covariate, and loci files for gemma runs and filters the dataset to make them.
    
    Parameters
    ----------
    splabs: list of strings
    Each string is the PATRIC taxon id for the species of one of the isolates. 
    Should be in same order as mutlist, phens; true if generated from other functions here. 
    mutlist: list of matrices
    Binary-encoded matrix for each gene of interest, one per gene of interest, with a row per isolate and col per site
    phens: list of ints
    0 or 1 depending on phenotype of interest, onne per isolate. 
    reflist: list of matrices
    Each matrix has the AA sequences for the reference isolates for each species, one per gene of interest
    genenames: list of strings
    Each string is the name of one of the genes of interest
    nametoid: dict
    Converts from text names of bacteria species in analysis (i.e. "M. tb") to id corresponding to 
    the associated row in reflist for that species
    anti: string
    Name of the antibiotic of interest
    species: string
    Name (i.e. "M. tb") of the species of interest. If None (default), all species with data will be included
    ref: string
    Name (i.e. "M. tb") of species whose reference blanks will be dropped from analysis. 
    Default is M. tb, overrides to be species indicated with species variable if not None. 
    
    Returns
    -------
    None, but makes three output files
    
    """
    
    ## FILTER BY SPECIES, IF APPLICABLE
    newmutlist = mutlist.copy()
    sp = nametoid[ref]
    sptype = "all" if species == None else species
    ## MAKE INDICATORS, MERGE MATRICES
    indlist, speciesdroplist, nonspidx = [], [], []
    indmat = []
    if not species: 
        for idx, spc in enumerate(["E. coli", "M. tb", "S. enterica"]):
            tempind = np.asarray([[1] if i == spc else [0] for i in splabs])
            if np.sum(tempind) < 2:
                speciesdroplist.append(idx)
            indlist.append(tempind)
        indmat = np.concatenate(indlist, 1)
        indmat = np.delete(indmat, speciesdroplist, 1)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
        splabs = [j for i, j in enumerate(splabs) if i not in nonspidx]
        ids = [j for i, j in enumerate(ids) if i not in nonspidx]
        np.savetxt(f"gemma_inputs_final/{sptype}_{anti}.covariates", indmat, delimiter="\t")
    
    if species: 
        # convert speies into index
        sp = nametoid[species]
        # find indices without species, drop them
        nonspidx = [i for i, j in enumerate(splabs) if j != species]
        for i in range(len(mutlist)): 
            newmutlist[i] = np.delete(mutlist[i], nonspidx, 0)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
        splabs = [j for i, j in enumerate(splabs) if i not in nonspidx]
        ids = [j for i, j in enumerate(ids) if i not in nonspidx]
    # save phenotype output for gemma
    with open(f"gemma_inputs_final/{sptype}_{anti}.phenotypes", 'w') as myfile:
        wr = csv.writer(myfile, delimiter ='\n')
        wr.writerow(phens)
    with open(f"gemma_inputs_final/{sptype}_{anti}_gids.txt", 'w') as myfile:
        wr = csv.writer(myfile, delimiter ='\n')
        wr.writerow(ids)
    
    ## DROP BLANKS IN REFERENCE SEQUENCE, LOW VARIATION SITES (<0.1%)
    # iterate over each gene
    site_list, totaldrop = [], []
    sitenum = np.arange(sum([m.shape[1] for m in newmutlist]))
    for j in range(len(newmutlist)):
        # find sites with low variation
        numvar = sum(newmutlist[j])
        # lowvarsites = np.where(numvar < 2)[0] # newmutlist[j].shape[0] / 200)[0]
        lowvarsites = [] # lowvarsites.tolist()
        # find blanks in reference sequence
        refblanks = []
        for i in range(newmutlist[j].shape[1]):
            if reflist[j].matrix[sp,i] == '-': 
                refblanks.append(i)
        # drop the appropriate loci
        droploci = list(set(lowvarsites + refblanks))
        site_list.append(np.array(["{}.{}".format(genenames[j], i + 1 - len([1 for a in refblanks if a < i])) 
                  if i not in refblanks else "_" for i in range(newmutlist[j].shape[1])]))
        # site_drops = list(set(lowvarsites).difference(set(refblanks)))
        site_list[j] = np.delete(site_list[j], droploci)
        newmutlist[j] = np.delete(newmutlist[j], droploci, 1)
        totaldrop += [i + sum([m.shape[1] for m in mutlist[:j]]) for i in droploci]
    sitenum = np.delete(sitenum, totaldrop)
    muts = np.concatenate(newmutlist, 1)
    sites = np.concatenate(site_list)
            
    ## MAKE INDICATORS, MERGE MATRICES
    indlist = []
    if not species: 
        spclist = ["E. coli", "M. tb", "S. enterica"]
        if anti == "rifampicin": 
            spclist = ["E. coli"]
        for spc in spclist:
            indlist.append(np.asarray([[1] if i == spc else [0] for i in splabs]))
        inds = np.concatenate(indlist, 1)
        np.savetxt(f"gemma_inputs_final/{sptype}_{anti}.covariates", inds, delimiter="\t")
            
    ## FORMAT LOCI FILE
    major_allele_code='0'
    minor_allele_code='1'
    major_allele_string="\"A\""
    minor_allele_string="TRUE"
    fileobj = open(f"gemma_inputs_final/{sptype}_{anti}.loci", "w")
    for idx, position in np.ndenumerate(sites): ## USED TO BE SITENUM - ATTEMPTING MORE INFORMATIVE VERSION! 
        # slice corresponding to all strains
        vec = muts[:, idx]
        is_major_allele = np.equal(vec, 0)
        string = np.array([minor_allele_code] * len(vec), dtype=object)
        string[is_major_allele[:,0]] = major_allele_code
        fileobj.write(f'{position},{major_allele_string},{minor_allele_string},{",".join(string)}\n')


def fildrop_fp(splabs, mutlist, phens, reflist, genenames, nametoid, anti, ids, species = None, 
            ref = "M. tb"):
    """
    This function makes the phenotype, covariate, and loci files for gemma runs for false-positive analyses
    
    Parameters
    ----------
    splabs: list of strings
    Each string is the PATRIC taxon id for the species of one of the isolates. 
    Should be in same order as mutlist, phens; true if generated from other functions here. 
    mutlist: list of matrices
    Binary-encoded matrix for each gene of interest, one per gene of interest, with a row per isolate and col per site
    phens: list of ints
    0 or 1 depending on phenotype of interest, onne per isolate. 
    reflist: list of matrices
    Each matrix has the AA sequences for the reference isolates for each species, one per gene of interest
    genenames: list of strings
    Each string is the name of one of the genes of interest
    nametoid: dict
    Converts from text names of bacteria species in analysis (i.e. "M. tb") to id corresponding to 
    the associated row in reflist for that species
    anti: string
    Name of the antibiotic of interest
    species: string
    Name (i.e. "M. tb") of the species of interest. If None (default), all species with data will be included
    ref: string
    Name (i.e. "M. tb") of species whose reference blanks will be dropped from analysis. 
    Default is M. tb, overrides to be species indicated with species variable if not None. 
    
    Returns
    -------
    None, but makes three output files
    
    """
    
    ## FILTER BY SPECIES, IF APPLICABLE
    newmutlist = mutlist.copy()
    sp = nametoid[ref]
    sptype = "all" if species == None else species
    ## MAKE INDICATORS, MERGE MATRICES
    indlist, speciesdroplist, nonspidx = [], [], []
    indmat = []
    if not species: 
        for idx, spc in enumerate(["E. coli", "M. tb", "S. enterica"]):
            tempind = np.asarray([[1] if i == spc else [0] for i in splabs])
            if np.sum(tempind) < 2:
                speciesdroplist.append(idx)
            indlist.append(tempind)
        indmat = np.concatenate(indlist, 1)
        indmat = np.delete(indmat, speciesdroplist, 1)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
        splabs = [j for i, j in enumerate(splabs) if i not in nonspidx]
        ids = [j for i, j in enumerate(ids) if i not in nonspidx]
        np.savetxt(f"gemma_inputs_final/fp_{sptype}_{anti}.covariates", indmat, delimiter="\t")
        
    if species: 
        # convert speies into index
        sp = nametoid[species]
        # find indices without species, drop them
        nonspidx = [i for i, j in enumerate(splabs) if j != species]
        for i in range(len(mutlist)): 
            newmutlist[i] = np.delete(mutlist[i], nonspidx, 0)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
        splabs = [j for i, j in enumerate(splabs) if i not in nonspidx]
        ids = [j for i, j in enumerate(ids) if i not in nonspidx]
    # save phenotype output for gemma
    with open(f"gemma_inputs_final/fp_{sptype}_{anti}.phenotypes", 'w') as myfile:
        wr = csv.writer(myfile, delimiter ='\n')
        wr.writerow(phens)
    with open(f"gemma_inputs_final/fp_{sptype}_{anti}_gids.txt", 'w') as myfile:
        wr = csv.writer(myfile, delimiter ='\n')
        wr.writerow(ids)
    
    ## DROP BLANKS IN REFERENCE SEQUENCE, LOW VARIATION SITES (<0.1%)
    # iterate over each gene
    site_list, totaldrop = [], []
    sitenum = np.arange(sum([m.shape[1] for m in newmutlist]))
    for j in range(len(newmutlist)):
        # find sites with low variation
        numvar = sum(newmutlist[j])
        # lowvarsites = np.where(numvar < 1)[0] # newmutlist[j].shape[0] / 200)[0]
        lowvarsites = [] #lowvarsites.tolist()
        # find blanks in reference sequence
        refblanks = []
        for i in range(newmutlist[j].shape[1]):
            if reflist[j].matrix[sp,i] == '-': 
                refblanks.append(i)
        # drop the appropriate loci
        droploci = list(set(lowvarsites + refblanks))
        site_list.append(np.array(["{}.{}".format(genenames[j], i + 1 - len([1 for a in refblanks if a < i])) 
                  if i not in refblanks else "_" for i in range(newmutlist[j].shape[1])]))
        # site_drops = list(set(lowvarsites).difference(set(refblanks)))
        site_list[j] = np.delete(site_list[j], droploci)
        newmutlist[j] = np.delete(newmutlist[j], droploci, 1)
        totaldrop += [i + sum([m.shape[1] for m in mutlist[:j]]) for i in droploci]
    sitenum = np.delete(sitenum, totaldrop)
    muts = np.concatenate(newmutlist, 1)
    sites = np.concatenate(site_list)
            
    ## FORMAT LOCI FILE
    major_allele_code='0'
    minor_allele_code='1'
    major_allele_string="\"A\""
    minor_allele_string="TRUE"
    fileobj = open(f"gemma_inputs_final/fp_{sptype}_{anti}.loci", "w")
    for idx, position in np.ndenumerate(sites): ## USED TO BE SITENUM - ATTEMPTING MORE INFORMATIVE VERSION! 
        # slice corresponding to all strains
        vec = muts[:, idx]
        is_major_allele = np.equal(vec, 0)
        string = np.array([minor_allele_code] * len(vec), dtype=object)
        string[is_major_allele[:,0]] = major_allele_code
        fileobj.write(f'{position},{major_allele_string},{minor_allele_string},{",".join(string)}\n')

## generates table with information about all sites in analysis
def table_info(mutlist, sites, sitenums, anti, genes, species = None, splabs = None): 
    muts = np.concatenate(mutlist, 1)
    counts = []
    if species: 
        ind = np.asarray([1 if i == species else 0 for i in splabs])
        counts = list(np.sum(muts[ind == 1], axis = 0))
    else: 
        counts = list(np.sum(muts, axis = 0))
    df = pd.DataFrame(list(zip(sites, sitenums, counts)), 
               columns =['Site', 'SiteNum', f"num {species if species != None else 'all'}"])
    df.to_csv(f"gemma_final_outputs/{anti}_{species if species != None else 'all'}_siteinfo.csv")

## generates table with information about all sites in false positive analysis
def table_info_fp(mutlist, sites, sitenums, anti, genes, species = None, splabs = None): 
    muts = np.concatenate(mutlist, 1)
    counts = []
    if species: 
        ind = np.asarray([1 if i == species else 0 for i in splabs])
        counts = list(np.sum(muts[ind == 1], axis = 0))
    else: 
        counts = list(np.sum(muts, axis = 0))
    df = pd.DataFrame(list(zip(sites, sitenums, counts)), 
               columns =['Site', 'SiteNum', f"num {species if species != None else 'all'}"])
    df.to_csv(f"gemma_final_outputs/fp_{anti}_{species if species != None else 'all'}_siteinfo.csv")

## calls the above functions for a given set of analyses
def pre_grm_gemma_prep(goi_alns, gois, gids, anti, nametoid, species_list, ref = "M. tb"): 
    # call aln_combiner to get 
    mutlist, phens, inds, reflist, idlist = aln_combiner(goi_alns, gids, anti)
    for species in species_list: 
        muts, sitenum, sites = fildrop(inds, mutlist, phens, reflist, gois, nametoid, anti, idlist, species, ref)
        table_info(mutlist, sites, sitenum, anti, gois, species, inds)

## calls the above functions for a given set of false positive analyses
def pre_grm_gemma_prep_fp(goi_alns, gois, gids, anti, nametoid, species_list, ref = "M. tb"): 
    # call aln_combiner to get 
    mutlist, phens, inds, reflist, idlist = aln_combiner(goi_alns, gids, anti)
    for species in species_list: 
        muts, sitenum, sites = fildrop_fp(inds, mutlist, phens, reflist, gois, nametoid, anti, idlist, species, ref)
        table_info_fp(mutlist, sites, sitenum, anti, gois, species, inds)

# set of dictionaries with conversions btwn different species referal conventions
sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}
sptoname = {"562": "E. coli", "1773": "M. tb", "1733": "M. tb", "28901": "S. enterica", "1280": "S. aureus"}
nametoid = {"S. aureus": 3, "E. coli": 0, "M. tb": 1, "S. enterica": 2}
gids = "final_combined_genomes.csv"

### RUN-SPECIFICS (CALLING FUNCTIONS)
## all species analysis, cipro
## specify variables
# list of strings of file name / pathways for alignments of genes of interest
goi_alns = ["afas/gyra.afa", "afas/gyrb.afa"]
# list of names of genes of interest (strings)
gois = ["gyrA", "gyrB"]
# string with name of antibiotic
anti = "ciprofloxacin"
# a list of strings with genome alignment file paths
species_list = [None, "E. coli", "S. enterica"]
# string with reference default
ref = "E. coli"
## call the function
pre_grm_gemma_prep(goi_alns, gois, gids, anti, nametoid, species_list, ref)

### RUN-SPECIFICS (CALLING FUNCTIONS)
## all species analysis, cip FALSE POSITIVE
## specify variables
# list of strings of file name / pathways for alignments of genes of interest
goi_alns = ["afas/rpoc.afa"]
# list of names of genes of interest (strings)
gois = ["rpoC"]
## call the function
pre_grm_gemma_prep_fp(goi_alns, gois, gids, anti, nametoid, species_list, ref)

### RUN-SPECIFICS (CALLING FUNCTIONS)
## all species analysis, rif
## specify variables
# list of strings of file name / pathways for alignments of genes of interest
goi_alns = ["afas/rpoa.afa", "afas/rpob.afa", "afas/rpoc.afa"]
# list of names of genes of interest (strings)
gois = ["rpoA", "rpoB", "rpoC"]
# string with name of antibiotic
anti = "rifampicin"
# a list of strings with genome alignment file paths
species_list = [None, "M. tb", "S. aureus"]
# string with reference default
ref = "M. tb"
## call function
pre_grm_gemma_prep(goi_alns, gois, gids, anti, nametoid, species_list, ref)

### RUN-SPECIFICS (CALLING FUNCTIONS)
## all species analysis, rif FALSE POSITIVE
## specify variables
# list of strings of file name / pathways for alignments of genes of interest
goi_alns = ["afas/gyrb.afa"]
# list of names of genes of interest (strings)
gois = ["gyrB"]
## call function
pre_grm_gemma_prep_fp(goi_alns, gois, gids, anti, nametoid, species_list, ref)
