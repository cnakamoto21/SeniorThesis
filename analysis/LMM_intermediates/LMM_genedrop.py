## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import csv

## SOME UNIVERSAL STUFF!
# make a dictionary encoding the order of the species
sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}
sptoname = {"562": "E.coli", "1773": "M.tb", "1733": "M.tb", "28901": "S.enterica", "1280": "S.aureus"}
speciestosp = {"S aureus": 3, "E.coli": 0, "M.tb": 1, "S.enterica": 2}

## reads in alignments and phenotype data and produces list of phenotypes based on alignment and antibioitc phenotype
def phen_cleaner(alignmentfile, rsstdata, antibiotic): 

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

## generrates single binary matrix from alignment
def aln_cleaner(alignmentfile, keepers = None, oldaln = None):
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
    # get reference sequences as an alignment object
    refs = aln.select(sequences = [aln.id_to_index["562.ref"], aln.id_to_index["1773.ref"], 
                                   aln.id_to_index["28901.ref"], aln.id_to_index["1280.ref"]])
    # create a matrix to conduct the identity screen
    identities = np.vstack((alncln.identities_to(refs[0]), alncln.identities_to(refs[1]), 
                            alncln.identities_to(refs[2]), alncln.identities_to(refs[3]))).T
    # make a new matrix of major and minor allele binaries
    muts = np.zeros((alncln.matrix.shape[0], alncln.matrix.shape[1]))
    muts.fill(np.nan)
    # iterate over the rows, initiate binary matrix and check identity level
    identityfails = []
    for i in range(muts.shape[0]):
        muts[i,:] = (refs.matrix[sptoid[specieslabs[i]]] != alncln.matrix[i,:]) * 1
        if identities[i,sptoid[specieslabs[i]]] < 0.9:
            identityfails.append(i)
    
    ## RETURN RELEVANT DATA
    return muts, alncln, refs, identityfails

## calls phenotype and alignment cleaners to generate binary matrices, phenotypes, labels for all genes
def aln_combiner(alignmentfiles, rsstdata, antibiotic):
    ## RUN OTHER CUSTOM FUNCTIONS TO GET THE ALIGNMENTS, PHENOTYPES
    # get phenotype data
    phens, keepers = phen_cleaner(alignmentfiles[0], rsstdata, antibiotic)
    mutlist, alnlist, reflist, fullidlist = [], [], [], []
    # iterate over the genes, run the alignment cleaner for each
    for i in range(len(alignmentfiles)):
        if i == 0:
            muts, alncln, refs, idfails = aln_cleaner(alignmentfiles[i], keepers)
            fullidlist = alncln.ids
        else: 
            muts, alncln, refs, idfails = aln_cleaner(alignmentfiles[i], oldaln = alnlist[0])
        # drop appropriate sites and add them to lists
        mutlist.append(np.delete(muts, idfails, 0))
        alnlist.append(alncln.select(sequences = [i for i in range(alncln.N) if i not in idfails]))
        reflist.append(refs)
    
    ## RESTRICT ISOLATES IN MUTS, PHENOTYPES
    # get a list of ids
    idsinall = alnlist[0].ids
    for aln in alnlist:
        idsinall = list(set(idsinall) & set(aln.ids))
    # find the ids not in the final list, drop at corresponding indices from matrices and alignments
    for i in range(len(mutlist)):
        drops = [alnlist[i].id_to_index[x] for x in alnlist[i].ids if x not in idsinall]
        mutlist[i] = np.delete(mutlist[i], drops, 0)
    # fix phenotype list
    drops = [i for i, x in enumerate(fullidlist) if x not in idsinall]
    phens = [j for i, j in enumerate(phens) if i not in drops]
    
    # CREATE SPECIES LABELS
    idlist = alnlist[0].select(sequences = [alnlist[0].id_to_index[x] for x in alnlist[0].ids if x in idsinall]).ids
    splabs = [sptoname[i.split(".")[0]] for i in idlist]
    
    return mutlist, phens, splabs, reflist

## CONSTRUCT THE THREE DIFFERENT TYPES OF GRMS!!
def grmer(splabs, mutlist, phens, reflist, genenames, anti, species = None, sp = 1, refsp = 1):

    ## FILTER BY SPECIES, IF APPLICABLE
    newmutlist = mutlist.copy()
    newmutlist = [m.astype(np.int8) for m in newmutlist]
    if species: 
        sp = speciestosp[species]
        nonspidx = [i for i, j in enumerate(splabs) if j != species]
        for i in range(len(mutlist)): 
            newmutlist[i] = np.delete(mutlist[i], nonspidx, 0)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
        splabs = [j for i, j in enumerate(splabs) if i not in nonspidx]
    with open(f"gemma_trial/{anti}_genedrop_{species if species != None else 'all'}.phenotypes", 'w') as myfile:
        wr = csv.writer(myfile, delimiter ='\n')
        wr.writerow(phens)
    
    ## DROP BLANKS IN REFERENCE SEQUENCE, LOW VARIATION SITES (<0.1%)
    # iterate over each gene
    sites, totaldrop = [], []
    sitenum = np.arange(sum([m.shape[1] for m in newmutlist]))
    for j in range(len(newmutlist)):
        # find sites with low variation
        numvar = sum(newmutlist[j])
        # lowvarsites = np.where(numvar < newmutlist[j].shape[0] / 1000)[0]
        lowvarsites = [] # lowvarsites.tolist()
        # find blanks in reference sequence
        refblanks = []
        for i in range(newmutlist[j].shape[1]):
            if reflist[j].matrix[sp,i] == '-': 
                refblanks.append(i)
        siterefblanks = refblanks
        if sp != refsp: 
            siterefblanks = []
            for i in range(newmutlist[j].shape[1]):
                if reflist[j].matrix[refsp,i] == '-': 
                    siterefblanks.append(i)
        # drop the appropriate loci
        droploci = list(set(lowvarsites + refblanks))
        # sites += ["{}.{}".format(genenames[j], i - len([1 for a in refblanks if a < i])) 
        #           for i in range(newmutlist[j].shape[1]) if i not in refblanks]
        site_temp = np.array(["{}.{}".format(genenames[j], i + 1 - len([1 for a in siterefblanks if a < i])) 
                  if i not in siterefblanks else "_" for i in range(newmutlist[j].shape[1])])
        site_temp = np.delete(site_temp, droploci)
        sites += site_temp.tolist()
        newmutlist[j] = np.delete(newmutlist[j], droploci, 1)
        totaldrop += [i + sum([m.shape[1] for m in mutlist[:j]]) for i in droploci]
    sitenum = np.delete(sitenum, totaldrop)
    muts = np.concatenate(newmutlist, 1)
                                     
    ## BUILD THE GRM - drop the gene of interest!!
    if species: 
        tempmuts = np.concatenate(newmutlist[1:], 1)
        # p = np.mean(tempmuts, axis = 0)
        # tempmuts = tempmuts - np.tile(p, (tempmuts.shape[0],1))
        grm = np.matmul(tempmuts, np.transpose(tempmuts)) / tempmuts.shape[1] # / sum(2 * p * (np.ones(p.shape) - p))
    else: 
        grm = np.zeros((mutlist[0].shape[0], mutlist[0].shape[0]))
        # iterate over each species
        for spcs in ["E.coli", "M.tb", "S.enterica", "S.aureus"]:
            # create a matrix with only the species
            temps = mutlist.copy()
            nonspidx = [i for i, j in enumerate(splabs) if j != spcs]
            for j in range(len(mutlist)): 
                temps[j] = np.delete(mutlist[j], nonspidx, 0)
            tempmuts = np.concatenate(temps[1:], 1)
            # p = np.mean(tempmuts, axis = 0)
            # tempmuts = tempmuts - np.tile(p, (tempmuts.shape[0],1))
            # get the pcs for the grm of the species
            minigrm = np.matmul(tempmuts, np.transpose(tempmuts)) / tempmuts.shape[1] # / sum(2 * p * (np.ones(p.shape) - p))
            # add that data into the overall matrix
            acnt = 0
            for a in range(mutlist[0].shape[0]): 
                if splabs[a] == spcs:
                    bcnt = 0
                    for b in range(mutlist[0].shape[0]):
                        if splabs[b] == spcs:
                            grm[a,b] = minigrm[acnt, bcnt]
                            bcnt += 1
                    acnt += 1
    np.savetxt(f"gemma_trial/{anti}_genedrop_{species if species != None else 'all'}.grm.txt", grm, delimiter="\t")
    
    ## MAKE INDICATORS, MERGE MATRICES, SAVE
    indlist = [np.ones((len(splabs), 1))]
    if not species:
        for spc in ["E.coli", "M.tb", "S.enterica"]:
            tempind = np.asarray([[1] if i == spc else [0] for i in splabs])
            if sum(tempind > 0):
                indlist.append(tempind)
        inds = np.concatenate(indlist, 1)
        np.savetxt(f"gemma_trial/{anti}_genedrop_{species if species != None else 'all'}.covariates", inds, delimiter="\t")
            
    ## FORMAT LOCI FILE
    major_allele_code='0'
    minor_allele_code='1'
    major_allele_string="\"A\""
    minor_allele_string="TRUE"
    fileobj = open(f"gemma_trial/{anti}_genedrop_{species if species != None else 'all'}.loci", "w")
    for idx, position in np.ndenumerate(sites):
        # slice corresponding to all strains
        vec = muts[:, idx]
        is_major_allele = np.equal(vec, 0)
        string = np.array([minor_allele_code] * len(vec), dtype=object)
        string[is_major_allele[:,0]] = major_allele_code
        fileobj.write(f'{position},{major_allele_string},{minor_allele_string},{",".join(string)}\n')
                      
    ## RETURN RELEVANT DATA
    return muts, sitenum, sites

## generate csv of relevant information about each site in the dataset for later analysis
def table_info(mutlist, sites, sitenums, anti, genes, species = None, splabs = None): 
    counts = []
    if species: 
        ind = np.asarray([1 if i == species else 0 for i in splabs])
    muts = np.concatenate(mutlist, 1)
    for idx in sitenums:
        if species: 
            counts.append(sum(muts[ind == 1][:,idx]))
        else: 
            counts.append(sum(muts[:,idx]))
    df = pd.DataFrame(list(zip(sites, sitenums, counts)), 
               columns =['Site', 'SiteNum', f"num {species if species != None else 'all'}"])
    df.to_csv(f"gemma_trial/{anti}_dropsite_{species if species != None else 'all'}_siteinfo.csv")

# individual species analysis, cipro
gyralnfiles = ["afas/gyra.afa", "afas/gyrb.afa"]
gyrmutlist, gyrphens, gyrinds, gyrreflist = aln_combiner(gyralnfiles, "combined_genomes.csv", "ciprofloxacin")
gyrgenes = ["gyrA", "gyrB"]
allgyrmuts, allgyrsitenums, allgyrsites = grmer(gyrinds, gyrmutlist, gyrphens, gyrreflist, gyrgenes,
                                                 "cip", species = "E.coli", sp = 0, refsp = 0)
table_info(gyrmutlist, allgyrsites, allgyrsitenums, "cip", gyrgenes, species = "E.coli", splabs = gyrinds)
allgyrmuts, allgyrsitenums, allgyrsites = grmer(gyrinds, gyrmutlist, gyrphens, gyrreflist, gyrgenes,
                                                 "cip", species = "S.enterica", sp = 0, refsp = 0)
table_info(gyrmutlist, allgyrsites, allgyrsitenums, "cip", gyrgenes, species = "S.enterica", splabs = gyrinds)

# individual species analysis, rif
rifalnfiles = ["afas/rpob.afa", "afas/rpoa.afa", "afas/rpoc.afa"]
rifmutlist, rifphens, rifinds, rifreflist = aln_combiner(rifalnfiles, "combined_genomes.csv", "rifampicin")
# note: for new code, must always include the gene of interest first. 
rifgenes = ["rpoB", "rpoA", "rpoC"]
allrifmuts, allrifsitenums, allrifsites = grmer(rifinds, rifmutlist, rifphens, rifreflist, rifgenes,
                                                 "rif", species = "M.tb", sp = 1, refsp = 1)
table_info(rifmutlist, allrifsites, allrifsitenums, "rif", rifgenes, species = "M.tb", splabs = rifinds)
allrifmuts, allrifsitenums, allrifsites = grmer(rifinds, rifmutlist, rifphens, rifreflist, rifgenes,
                                                 "rif", species = "S.aureus", sp = 3, refsp = 1)
table_info(rifmutlist, allrifsites, allrifsitenums, "rif", rifgenes, species = "S.aureus", splabs = rifinds)
