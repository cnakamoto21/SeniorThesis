## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

## SOME UNIVERSAL STUFF!
# make a dictionary encoding the order of the species
sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}
speciestoid = {"E. coli": 0, "M. tb": 1, "S. enterica": 2, "S. aureus": 3}

## phenotype cleaner: reads in phenotypes, generates list in same order as alignment
## also identifies isolates without phenotype of interest
def phen_cleaner(alignmentfile, rsstdata, antibiotic): 
    
    ## READ IN RELEVANT DATA
    # get an alignment object (with matrix and stuff) from alignment file
    with open(alignmentfile, "r") as infile:
        aln = Alignment.from_file(infile, format="fasta")
    # get all phenotype data
    genomepd = pd.read_csv('combined_genomes.csv', dtype = str)

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

def aln_cleaner(alignmentfile, keepers = None, oldaln = None):
    """
    This function sequence alignments for a single gene w/ resistance data for each isolate. 
    1) reads in the alignment file and resistance data. 
    2) recodes ressitance data as binary and reorders to fit alignment order. 
    3) restricts alignment data based on phenotype data availability. 
    4) recodes alignment 1 at positions that deviate from species reference else 0.
    5) one hot encodes the species labels and adds them to the alignment data. 
    6) applies 90% identity filter to the amino acid sequences wrt species reference
    Note: function assumes data from E. coli, M. tuberculosis, S. aureus, and S. enterica
    
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
    2D numpy array
        a row for each isolate and a column for each aa position. cleaned binary matrix (see 4-6 above).
    list of integers
        indicates susceptible (0) or resistant (1) phenotypes in same order as numpy array rows. 
    ev couplings Alignment
        the raw ev couplings alignment file corresponding to the file
    list of integers
        the indices for the e. coli, tb, salmonella, and staph reference sequences (in order) in the Alignment
    """
    
    ## IMPORT RELEVANT PACKAGES / FUNCTIONS
    from evcouplings.align import Alignment, map_matrix
    import numpy as np
    
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

## runs alignment and phenotype cleaner functions to generate cleaned binary matrices, phenotype list,
## array of indicators, and references for each species
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
    
    # CREATE INDICATORS
    idlist = alnlist[0].select(sequences = [alnlist[0].id_to_index[x] for x in alnlist[0].ids if x in idsinall]).ids
    specieslabs = [i.split(".")[0] for i in idlist]
    ecoliind = np.asarray([[1] if i == '562' else [0] for i in specieslabs])
    tbind = np.asarray([[1] if i == '1773' or i == '1733' else [0] for i in specieslabs])
    salmind = np.asarray([[1] if i == '28901' else [0] for i in specieslabs])
    indicators = np.concatenate((ecoliind, tbind, salmind), 1)
    
    return mutlist, phens, indicators, reflist

def lg_reg_grm(indicators, mutlist, phens, reflist, genenames, species = None, sp = 1, refsp = 1, npca = 10):
    """
    This function takes the outputs of aln_cleaner and performs logistic regressions for each site.
    It also first filters the variant matrix of some sites and isolates out of the indicated species if applicable.
    It also generates a GRM from the variant matrix and does the PCA transformation of the GRM, incorporating PCs into the logistic regressions. 
    
    Parameters
    ----------
    indicators: 2D numpy matrix
    The thrid output of aln_cleaner, indicators for the species of each isolate
    mutlist: list of 2D numpy matrices
    The first output from aln_cleaner, list of binary encoded variant matrices for each gene
    phens: list of ints (binary)
    The second output from aln_cleaner with the list of phenotypes, 1 for resistant and 0 for susceptible
    reflist: list of integers
    The fourth output from aln_cleaner, list of the indices for the reference sequences in the Alignment
    genenames: list of strings
    List of the names of all of the genes included in the analysis, in the same order in which their alignment files were fed into aln_cleaner
    species: string
    Name of the species. Options are "E. coli", "M. tuberculosis", "S. enterica", "S. aureus", and None for all-species
    sp: int
    Indicates the reference species in the case of an all-species analysis else overriden by the species variable. 
    See accompanying dicts for int options (0 for E. coli, 1 for M. tuberculosis, etc.)
    refsp: int
    Indicates the species to use for site naming. See accompanying dicts for int options (0 for E. coli, 1 for M. tuberculosis, etc.)
    npca: int
    The number of principal components to include in the regression from each species in the analysis

    
    Returns
    -------
    list of ints
        the p-value for the amino acid position of interest in each logistic regression
    2D numpy array
        a row for each isolate and a column for each aa position. cleaned binary matrix.
    2D numpy array
        same as above but without low-variation sites or sites that are blanks in the reference. 
    1D numpy array
        the names of all of the sites in the amino acid sequence remaining after cleaning, in order
    1D numpy array
        the indices of all of the sites in the amino acid sequence remaining after cleaning in input muts
    list of ints
        the odds ratio for the amino acid position of interest in each logistic regression   
    """
    
    ## IMPORT RELEVANT PACKAGES / FUNCTIONS
    from evcouplings.align import Alignment, map_matrix
    import numpy as np
    import statsmodels.api as sm
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    ## FILTER BY SPECIES, IF APPLICABLE
    newmutlist = [np.copy(m).astype(np.int8) for m in mutlist]
    # newmutlist = mutlist.copy()
    # newmutlist = [m.astype(np.int8) for m in newmutlist]
    speciesdroplist = []
    if species: 
        if species == "S. aureus": 
            sp = 3
            indicator = np.ones(mutlist[0].shape[0]) - np.sum(indicators, 1)
        else:
            if species == "E. coli":
                sp = 0
            elif species == "M. tb":
                sp = 1
            else: 
                sp = 2
            indicator = indicators[:,sp]
        nonspidx = [i for i in range(mutlist[0].shape[0]) if indicator[i] == 0]
        for i in range(len(mutlist)): 
            newmutlist[i] = np.delete(newmutlist[i], nonspidx, 0)
        indicators = np.delete(indicators, nonspidx, 0)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx] 
        # get the 4 pcs of the grms
        tempmuts = np.concatenate(newmutlist, 1)
        return tempmuts, tempmuts, tempmuts, tempmuts, tempmuts
        grm = np.matmul(tempmuts, np.transpose(tempmuts)) / tempmuts.shape[1]
        pcs = PCA(n_components = npca).fit_transform(grm) 
    else: 
        pcs = np.zeros((newmutlist[0].shape[0], 4 * npca))
        # iterate over each species
        for i in range(4):
            if i == 3:
                ind = np.ones(mutlist[0].shape[0]) - np.sum(indicators, 1)
            else:
                ind = indicators[:,i]
            if np.sum(ind) < 10: 
                pcs = np.delete(pcs, list(range(i*npca,(i+1)*npca)), 1)
                if i < 3: 
                    speciesdroplist.append(i)
                continue
            # create a matrix with only the species
            temps = mutlist.copy()
            nonspidx = [i for i in range(mutlist[0].shape[0]) if ind[i] == 0]
            for j in range(len(mutlist)): 
                temps[j] = np.delete(mutlist[j], nonspidx, 0)
            tempmuts = np.concatenate(temps, 1)
            # get the pcs for the grm of the species
            grm = np.matmul(tempmuts, np.transpose(tempmuts)) / tempmuts.shape[1]
            pc4 = PCA(n_components = npca).fit_transform(grm)
            # add that data into the overall matrix
            pccounter = 0
            for a in range(newmutlist[0].shape[0]): 
                if ind[a] == 1:
                    pcs[a,i:i + npca] = pc4[pccounter]
                    pccounter += 1
        
    ## DROP BLANKS IN REFERENCE SEQUENCE, LOW VARIATION SITES (<0.1%)
    # iterate over each gene
    sites, totaldrop = [], []
    sitenum = np.arange(sum([m.shape[1] for m in newmutlist]))
    for j in range(len(newmutlist)):
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
        # find sites with low variation
        numvar = sum(newmutlist[j])
        lowvarsites = np.where(numvar < newmutlist[j].shape[0] / 1000)[0] # ADJUST THIS AS NECESSARY!!
        lowvarsites = lowvarsites.tolist()
        # drop the appropriate loci
        # droploci = list(set(lowvarsites + refblanks))
        # sites += ["{}.{}".format(genenames[j], i) for i in range(newmutlist[j].shape[1]) if i not in droploci]
        site_temp = np.array(["{}.{}".format(genenames[j], i + 1 - len([1 for a in siterefblanks if a < i])) 
                  if i not in siterefblanks else "_" for i in range(newmutlist[j].shape[1])])
        # site_temp = np.delete(site_temp, droploci)
        sites += site_temp.tolist()
        # newmutlist[j] = np.delete(newmutlist[j], droploci, 1)
        # totaldrop += [i + sum([m.shape[1] for m in mutlist[:j]]) for i in droploci]
    # sitenum = np.delete(sitenum, totaldrop)
            
    ## MERGE THE MATRICES
    muts = np.concatenate(newmutlist, 1)
    # make another low var cut if the PCs are weird
    numvar = np.sum(pcs, axis = 0)
    lowvarsites = np.where(numvar == 0)[0]
    lowvarsites = lowvarsites.tolist()
    pcs = np.delete(pcs, lowvarsites, 1)
    indicators = np.delete(indicators, speciesdroplist, 1)
    
    ## CONDUCT LOGISTIC REGRESSION - OLS using StatsModels
    # function-specific data prep: convert response data, add intercept constants
    arphens = np.asarray(phens)
    betas, pvals = [], []
    if not species: 
        pcs = np.concatenate((pcs, indicators), axis = 1)
    for i in range(muts.shape[1] - 1): 
        X_added = sm.add_constant(np.concatenate((muts[:,i].reshape(-1,1), pcs), axis = 1))
        # conduct regression
        logreg = sm.Logit(arphens, X_added)
        # save regression info (parameters, etc) in results_sm
        results_sm = logreg.fit(method='bfgs', disp = False) # try nm
        # idx = 0 if species != None else 1
        idx = 1
        # get the relevant model outputs
        betas.append(np.exp(results_sm.params[idx]))
        pvals.append(results_sm.pvalues[idx])
        
    ## RETURN RELEVANT DATA
    return pvals, muts, sites, sitenum, betas

def qq_plot(p_val_list, title, filename, xlim=[0,10]):
    p_val_list = [p for p in p_val_list if not np.isnan(p)]
    sorted_p_vals = sorted(p_val_list)
    sorted_p_vals = [i if i != 0 else 1.0E-50 for i in sorted_p_vals]
    size = len(sorted_p_vals)
    
    increment = 1/size
    
    expected_p_vals = np.arange(0+increment,1+increment,increment)
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    
    expected = -np.log10(expected_p_vals)
    actual = -np.log10(sorted_p_vals)
    
    ax.scatter(expected, actual)
    ax.plot(expected, expected, color="black")
    ax.axhline(-np.log10(0.05 / size), color = "red")
    ax.set_ylim(xlim)
    ax.set_xlabel("expected pvalues")
    ax.set_ylabel("actual pvalues")
    ax.set_title(f'{title}')
    sns.despine()
    plt.savefig(f'figures/{filename}.png', dpi=600);
    plt.show()

def framer(pvals, sites, sitenum, reflist, betas, group, sp = 1):
    badpsites = [i for i, p in enumerate(pvals) if np.isnan(p)]
    pvals = [p for i, p in enumerate(pvals) if i not in badpsites]
    betas = [p for i, p in enumerate(betas) if i not in badpsites]
    sitenum = [p for i, p in enumerate(list(sitenum)) if i not in badpsites]
    sites = [p for i, p in enumerate(sites) if i not in badpsites]
    cutoffnum = sum([1 if i < (0.05 / len(pvals)) else 0 for i in sorted(pvals)])
    impsites = np.argsort(pvals).tolist()
    impnames = [sites[i] for i in impsites[0:cutoffnum]]
    impnums = [sitenum[i] for i in impsites[0:cutoffnum] if sitenum[i]]
    refs = np.concatenate([reflist[i].matrix for i in range(len(reflist))], 1)
    aaofinterest = [refs[sp,i] if i < refs.shape[1] else 'ind' for i in impnums]
    ps = sorted(pvals)[0:cutoffnum]
    beta = list(np.array(betas)[np.argsort(pvals)])
    frame = pd.DataFrame(list(zip(impnames, impnums, aaofinterest, beta, ps)), 
                                columns = ['Site', 'SiteNum', 'RegAA', f'OR-{group}', f'p-{group}']) 
    return frame, cutoffnum, impsites

def corranalysis(muts, cutoffnum, impsites, imp, title, filename):
    ## CORRELATION ANALYSIS
    # note: uses variable definitions from all species set
    importantsites = np.empty([cutoffnum, muts.shape[0]])
    for i in range(cutoffnum):
        importantsites[i,:] = muts[:,impsites[i]]
    corrs = np.corrcoef(importantsites)

    import seaborn as sns
    import matplotlib.pyplot as plt

    # plot the heatmap
    fig = plt.figure(figsize=(11,10))
    ax = fig.gca
    sns.heatmap(corrs, cmap="Blues", xticklabels=imp, yticklabels=imp).set_title(f'{title}')
    plt.savefig(f'figures/{filename}.png', dpi=600);
    
    return corrs

def frame_cleaner(frame, sitenumlists, species, plist, betalist, anti):
    for index, row in frame.iterrows():
        for i in range(len(species)):
            if np.isnan(row[f'p-{species[i]}']):
                if row['SiteNum'] in sitenumlists[i]: 
                    frame.loc[index, f'p-{species[i]}'] = "{:.1E}".format(plist[i][sitenumlists[i].index(row['SiteNum'])])
                    frame.loc[index, f'OR-{species[i]}'] = "{:.2f}".format(betalist[i][sitenumlists[i].index(row['SiteNum'])])
            elif isinstance(row[f'p-{species[i]}'], float): 
                frame.loc[index, f'p-{species[i]}'] = "{:.1E}*".format(frame.loc[index, f'p-{species[i]}'])
                frame.loc[index, f'OR-{species[i]}'] = "{:.2f}".format(frame.loc[index, f'OR-{species[i]}'])
    collist = ['Site'] + [f'num-{s}' for s in species] + [f'OR-{s}' for s in species] + [f'p-{s}' for s in species]
    frame = frame[collist]
    frame.to_csv(f'PCR_{anti}.csv', index = False)
    return frame

## Rifampicin Analyses
## all species analysis
# read in and clean data for analysis
alignmentfiles = ["afas/rpob.afa", "afas/rpoa.afa", "afas/rpoc.afa"]
rpomutlist, rpophens, rpoinds, rporeflist = aln_combiner(alignmentfiles, "combined_genomes.csv", "rifampicin")
# further clean and run logistic regression on data
genes = ["rpoB", "rpoA", "rpoC"]
allrpop, allrpomuts, allrposites, allrpositenums, allrpobetas = lg_reg_grm(rpoinds, rpomutlist, rpophens, 
                                                    rporeflist, genes, species = None, sp = 1, refsp = 1, npca = 3)
# visualize the results
allrpo, allrpopcut, allimpsites = framer(allrpop, allrposites, allrpositenums, rporeflist, allrpobetas, "all", 1)
qq_plot(allrpop, 'Multi-Species PC', 'PC_rif_all_qq')
allcorrs = corranalysis(allrpomuts, allrpopcut, allrpo.SiteNum, allrpo.Site, 'Multi-Species PC', 'PC_rif_all_corr')

## tb analysis
tbrpop, tbrpomuts, tbrposites, tbrpositenums, tbbetas = lg_reg_grm(rpoinds, rpomutlist, rpophens, 
                                                    rporeflist, genes, species = "M. tb", npca = 1)
# visualize results
tbrpo, tbrpopcut, tbrpoimpsites = framer(tbrpop, tbrposites, tbrpositenums, rporeflist, tbbetas, "tb", 1)
newrpob = pd.merge(allrpob, tbrpob, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(tbrpop, 'M. tuberculosis PC', 'PC_rif_tb_qq')
tbcorrs = corranalysis(tbrpomuts, tbrpopcut, tbrpoimpsites, tbrpo.Site, 'M. tuberculosis PC', 'PC_rif_tb_corr')

## s. aureus analysis
# further clean and run logistic regression on data
strpop, strpomuts, strposites, strpositenums, stbetas = lg_reg_grm(rpoinds, rpomutlist, rpophens, 
                                                            rporeflist, genes, species = "S. aureus", npca = 1)
# visualize results
strpo, strpopcut, strpoimpsites = framer(strpop, strposites, strpositenums, rporeflist, stbetas, "staph", 3)
newrpo2 = pd.merge(newrpo, strpo, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(strpop, 'S. aureus PC', 'PC_rif_st_qq')
if stpcut > 1: 
    stcorrs = corranalysis(strpomuts, strpopcut, strpoimpsites, strpo.Site, 'S. aureus PC', 'PC_rif_st_corr')

## merge and clean tables
numall = []
numtb = []
numstaph = []
mutfinder = {'rpoB': 0, "rpoA": 1, "rpoC":2}
for i in range(len(newrpo2)): 
    name = [mutfinder[newrpo2.loc[i,"Site"].split(".")[0]]]
    num = allrpositenums[allrposites.index(newrpo2.loc[i,"Site"])]
    name.append(num - sum([rpomutlist[i].shape[1] for i in range(name[0])]))
    numall.append(sum(rpomutlist[name[0]][:,name[1]]))
    numtb.append(sum(rpomutlist[name[0]][rpoinds[:,1] == 1][:,name[1]]))
    numstaph.append(sum(rpomutlist[name[0]][(rpoinds[:,0]==0)&(rpoinds[:,1]==0)&(rpoinds[:,2]==0)][:,name[1]]))
frame = pd.DataFrame(list(zip(newrpo2.Site, numall, numtb, numstaph)), 
                                columns = ['Site', 'num-all', 'num-tb', 'num-staph']) 
newrpo3 = pd.merge(newrpo2, frame, how = 'outer', on = ['Site'])
newrpo3 = newrpo3.astype({'num-all': 'int32', 'num-tb': 'int32', 'num-staph': 'int32'}) 
newrpo3['SiteNum'] = [allrpositenums[allrposites.index(i)] for i in newrpo3['Site']]
sitenumlists = [list(strpositenums), list(tbrpositenums), list(allrpositenums)]
species = ["tb", "staph", "all"] # 
plist = [tbrpop, strpop, allrpop] # 
betalist = [tbbetas, stbetas, allrpobetas] # 
newrpo4 = frame_cleaner(newrpo3, sitenumlists, species, plist, betalist, "rif")


## Ciprofloxacin Analyses
## all species analysis
# read in and clean data for analysis
alignmentfiles = ["afas/gyra.afa", "afas/gyrb.afa"]
gyrmutlist, gyrphens, gyrinds, gyrreflist = aln_combiner(alignmentfiles, "combined_genomes.csv", "ciprofloxacin")
# further clean and run logistic regression on data
gyrgenes = ["gyrA", "gyrB"]
allgyrpgrm, allgyrmutsgrm, allgyrsitesgrm, allgyrsitenumsgrm, allgyrbetas = lg_reg_grm(gyrinds, gyrmutlist, gyrphens, 
                                                gyrreflist, gyrgenes, species = None, sp = 0, refsp = 0, npca = 4)
# visualize the results
allgyrgrm, allgyrpcutgrm, allgyrimpsitesgrm = framer(allgyrpgrm, allgyrsitesgrm, allgyrsitenumsgrm, gyrreflist, 
                                                  allgyrbetas, "all", 0)
qq_plot(allgyrpgrm, 'Multi-Species PC', 'PC_cip_all_qq')
allgyrcorrs = corranalysis(allgyrmutsgrm, allgyrpcutgrm, allgyrimpsitesgrm, allgyrgrm.Site, 
                           'Multi-Species PC', 'PC_cip_all_corr')
## e. coli analysis
# further clean and run logistic regression on data
ecgyrpgrm, ecgyrmutsgrm, ecgyrsitesgrm, ecgyrsitenumsgrm, ecbetas = lg_reg_grm(gyrinds, gyrmutlist, gyrphens, 
                                                        gyrreflist, gyrgenes, species = "E. coli", refsp = 0, npca = 2)
# visualize results
ecgyrgrm, ecgyrpcutgrm, ecgyrimpsitesgrm = framer(ecgyrpgrm, ecgyrsitesgrm, ecgyrsitenumsgrm, gyrreflist, 
                                                  ecbetas, "ecoli", 0)
newgyrgrm = pd.merge(allgyrgrm, ecgyrgrm, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(ecgyrpgrm, 'E. coli PC', 'PC_cip_ec_qq')
ecgyrcorrs = corranalysis(ecgyrmutsgrm, ecgyrpcutgrm, ecgyrimpsitesgrm, ecgyrgrm.Site, 
                          'E. coli PC', 'PC_cip_ec_corr')

## s. enterica analysis
# further clean and run logistic regression on data
sagyrpgrm, sagyrmutsgrm, sagyrsitesgrm, sagyrsitenumsgrm, sabetas = lg_reg_grm(gyrinds, gyrmutlist, gyrphens, 
                                                    gyrreflist, gyrgenes, species = "S. enterica", refsp = 0, npca = 1)
# visualize results
sagyrgrm, sagyrpcutgrm, sagyrimpsitesgrm = framer(sagyrpgrm, sagyrsitesgrm, sagyrsitenumsgrm, gyrreflist, 
                                                  sabetas, "sal", 2)
newgyr2grm = pd.merge(newgyrgrm, sagyrgrm, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(sagyrpgrm, 'S. enterica PC', 'PC_cip_sal_qq')
if stpcut > 1: 
    sagyrcorrs = corranalysis(sagyrmutsgrm, sagyrpcutgrm, sagyrimpsitesgrm, sagyrgrm.Site, 'S. enterica PC', 'PC_cip_sal_corr')

## merge and clean tables
# get counts of mutants at each site
numall = []
numtb = []
numstaph = []
mutfinder = {'gyrA': 0, "gyrB": 1}
for i in range(len(newgyr2grm)): 
    if newgyr2grm.loc[i,"RegAA"] == "ind": 
        numall.append(0)
        numtb.append(0)
        numstaph.append(0)
        continue
    name = [mutfinder[newgyr2grm.loc[i,"Site"].split(".")[0]]]
    name.append(newgyr2grm.loc[i,"SiteNum"] - sum([gyrmutlist[i].shape[1] for i in range(name[0])]))
    numall.append(sum(gyrmutlist[name[0]][:,name[1]]))
    numtb.append(sum(gyrmutlist[name[0]][gyrinds[:,0] == 1][:,name[1]]))
    numstaph.append(sum(gyrmutlist[name[0]][gyrinds[:,2] == 1][:,name[1]]))
# add counts to dataframes, change type
frame = pd.DataFrame(list(zip(newgyr2grm.SiteNum, numall, numtb, numstaph)), 
                                columns = ['SiteNum', 'num-all', 'num-ecoli', 'num-sal']) 
newgyr3 = pd.merge(newgyr2grm, frame, how = 'outer', on = ['SiteNum'])
newgyr3 = newgyr3.astype({'num-all': 'int32', 'num-ecoli': 'int32', 'num-sal': 'int32'})
# fill in remaining blanks in merged table, where appropriate
## ADD THE OTHER OR, P-VALS BACK IN; FORMAT TABLE DISPLAY; SAVE
sitenumlists = [list(ecgyrsitenumsgrm), list(sagyrsitenumsgrm), list(allgyrsitenumsgrm)]
species = ['ecoli', 'sal', 'all']
plist = [ecgyrpgrm, sagyrpgrm, allgyrpgrm]
betalist = [ecbetas, sabetas, allgyrbetas]
anti = 'cip'
newgyr4 = frame_cleaner(newgyr3, sitenumlists, species, plist, betalist, anti)
