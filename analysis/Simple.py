## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt

def aln_cleaner(alignmentfile, rsstdata, antibiotic):
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
    import pandas as pd
    
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
    # drop the bad indices
    alncln = aln.select(sequences = keepers.astype(bool))
    
    ## CREATE BINARY MATRIX W/ SPECIES INDICATORS AND 90% SCREEN
    # make arrays with the indicators
    specieslabs = [i.split('.')[0] for i in list(alncln.ids)]
    # strain out incredibly underrepresented species
    for spcs in list(set(specieslabs)): 
        if sum([i == spcs for i in specieslabs]) < 10: 
            minispcs = [i != spcs for i in specieslabs]
            alncln = alncln.select(sequences = minispcs)
            for i in range(len(specieslabs)): 
                if specieslabs[i] == spcs: 
                    phens.pop(i)
    specieslabs = [i.split('.')[0] for i in list(alncln.ids)]
    # make indicators
    ecoliind = np.asarray([[1] if i == '562' else [0] for i in specieslabs])
    tbind = np.asarray([[1] if i == '1773' or i == '1733' else [0] for i in specieslabs])
    salmind = np.asarray([[1] if i == '28901' else [0] for i in specieslabs])
    # get reference sequence locations as list (order: e coli, tb, salmonella, staph)
    reflocs = [aln.id_to_index["562.ref"], aln.id_to_index["1773.ref"], 
               aln.id_to_index["28901.ref"], aln.id_to_index["1280.ref"]]
    # create a matrix to conduct the identity screen
    identities = np.vstack((alncln.identities_to(aln.matrix[reflocs[0]]), 
                            alncln.identities_to(aln.matrix[reflocs[1]]), 
                            alncln.identities_to(aln.matrix[reflocs[2]]), 
                            alncln.identities_to(aln.matrix[reflocs[3]]))).T
    # make a new matrix of major and minor alleles
    muts = np.zeros((alncln.matrix.shape[0], alncln.matrix.shape[1]))
    muts.fill(np.nan)
    # iterate over the rows, initiate binary matrix and check identity level
    identityfails = []
    for i in range(muts.shape[0]):
        if ecoliind[i,0] == 1:
            muts[i,:] = (aln.matrix[reflocs[0],:] != alncln.matrix[i,:]) * 1
            if identities[i,0] < 0.9:
                identityfails.append(i)
        elif tbind[i,0] == 1:
            muts[i,:] = (aln.matrix[reflocs[1],:] != alncln.matrix[i,:]) * 1
            if identities[i,1] < 0.9:
                identityfails.append(i)
        elif salmind[i,0] == 1:
            muts[i,:] = (aln.matrix[reflocs[2],:] != alncln.matrix[i,:]) * 1
            if identities[i,2] < 0.9:
                identityfails.append(i)
        else: 
            muts[i,:] = (aln.matrix[reflocs[3],:] != alncln.matrix[i,:]) * 1
            if identities[i,3] < 0.9:
                identityfails.append(i)
    # add the indicators to the matrix
    for col in [ecoliind, tbind, salmind]:
        muts = np.append(muts, col, 1)
    # cleanse relevant matrices of bad indices
    muts = np.delete(muts, identityfails, 0)
    phens = [j for i, j in enumerate(phens) if i not in identityfails]
    
    ## RETURN RELEVANT DATA
    return muts, phens, aln, reflocs

def lg_regressor(aln, muts, phens, reflocs, species = None, sp = 1, refsp = 1):
    """
    This function takes the outputs of aln_cleaner and performs logistic regressions for each site. 
    
    Parameters
    ----------
    aln: ev couplings Alignment
    The uncleaned ev couplings Alignment file for the gene of interest
    muts: 2D numpy matrix
    The first output from aln_cleaner
    phens: list of ints (binary)
    The second output from aln_cleaner with the list of phenotypes, 1 for resistant and 0 for susceptible
    reflocs: list of integers
    The fourth output from aln_cleaner, list of the indices for the reference sequences in the Alignment
    species: string
    Name of the species. Options are "E. coli", "M. tuberculosis", "S. enterica", "S. aureus", and None for all-species
    sp: int
    Indicates the reference species in the case of an all-species analysis else overriden by the species variable. 
    See accompanying dicts for int options (0 for e coli, 1 for tb, 2 for salmonella, 3 for staph)
    refsp: int
    Indicates the species to use for site naming. See accompanying dicts for int options (0 for E. coli, 1 for M. tuberculosis, etc.)
    
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
    
    ## FILTER BY SPECIES, IF APPLICABLE
    if species: 
        if species == "E. coli":
            sp = 0
            indicator = muts[:,-3]
        elif species == "M. tuberculosis":
            sp = 1
            indicator = muts[:,-2]
        elif species == "S. enterica": 
            sp = 2
            indicator = muts[:,-1]
        else: 
            sp = 3
            indicator = np.ones(muts.shape[0]) - np.sum(muts[:,-3:], 1)
        nonspidx = [i for i in range(muts.shape[0]) if indicator[i] == 0]
        muts = np.delete(muts, nonspidx, 0)
        phens = [j for i, j in enumerate(phens) if i not in nonspidx]
    
    ## DROP BLANKS IN REFERENCE SEQUENCE, LOW VARIATION SITES (<0.1%)
    indicators = muts[:,-3:]
    # find blanks in reference sequence
    refblanks = []
    for i in range(aln.matrix.shape[1]):
        if aln.matrix[reflocs[sp],i] == '-': 
            refblanks.append(i)
    # make a different reference set if necessary for site naming
    siterefblanks = refblanks
    if sp != refsp: 
        siterefblanks = []
        for i in range(aln.matrix.shape[1]):
            if aln.matrix[reflocs[refsp],i] == '-': 
                siterefblanks.append(i)
    # find sites with low variation
    numvar = sum(muts)
    lowvarsites = np.where(numvar < muts.shape[0] / 1000)[0]
    lowvarsites = lowvarsites.tolist()
    # drop the appropriate loci
    droploci = list(set(lowvarsites + refblanks))
    specmuts = np.delete(muts, droploci, 1)
    sites = np.array(["{}".format(i + 1 - len([1 for a in siterefblanks if a < i])) 
                  if i not in siterefblanks else "_" for i in range(aln.matrix.shape[1])])
    sites = np.delete(sites, droploci)
    # sites = sites.astype(np.int)
    # get the site number (for indexing into raw mutant matrix)
    sitenum = np.arange(aln.matrix.shape[1])
    sitenum = np.delete(sitenum, droploci)
    
    ## CONDUCT LOGISTIC REGRESSION - OLS using StatsModels
    # function-specific data prep: convert response data, add intercept constants
    arphens = np.asarray(phens)
    # drop unnecessary indicators first
    if species == None: 
        numvar = sum(indicators)
        badindicators = np.where(numvar == 0)[0]
        lowvarsites = badindicators.tolist()
        indicators = np.delete(indicators, lowvarsites, 1)
    betas, pvals = [], []
    for i in range(specmuts.shape[1] - 1): 
        if species: 
            X_added = sm.add_constant(specmuts[:,i])
        else: 
            X_added = sm.add_constant(np.concatenate((specmuts[:,i].reshape(-1,1), indicators), axis = 1))
        # conduct regression
        logreg = sm.Logit(arphens, X_added)
        # save regression info (parameters, etc) in results_sm
        results_sm = logreg.fit(method='bfgs', disp = False)
        # idx = 0 if species != None else 1
        idx = 1
        # get the relevant model outputs
        betas.append(np.exp(results_sm.params[idx]))
        pvals.append(results_sm.pvalues[idx])
    
    ## RETURN RELEVANT DATA
    return pvals, muts, specmuts, sites, sitenum, betas

def qq_plot(p_val_list, filename, title, xlim=[0,10]):
    ## generates png file qq plot from p-value data
    sorted_p_vals = sorted(p_val_list)
    size = len(sorted_p_vals)
    increment = 1/size
    expected_p_vals = np.arange(0+increment,1+increment,increment)
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    sorted_p_vals = [i if i != 0 else 1.0E-70 for i in sorted_p_vals]
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
    plt.savefig(f'{filename}.png', dpi=600);
    plt.show()

def framer(pvals, sites, sitenums, aln, reflocs, betas, sp, muts, group, extimp = []):
    ## generates dataframe with information about significant hits
    cutoffnum = sum([1 if i < (0.05 / len(pvals)) else 0 for i in sorted(pvals)])
    impsites = np.argsort(pvals).tolist()
    imp = sites[impsites[0:cutoffnum]].tolist()
    impsitenums = sitenums[impsites[0:cutoffnum]].tolist()
    aaofinterest = [aln.matrix[reflocs[sp],i] for i in impsitenums]
    ps = sorted(pvals)[0:cutoffnum]
    beta = list(np.array(betas)[np.argsort(pvals)])
    frame = pd.DataFrame(list(zip(imp, impsitenums, aaofinterest, beta, ps)), 
                                columns = ['Site', 'SiteNum', 'RegAA', f'OR-{group}', f'p-{group}']) 
    return frame, cutoffnum, impsites

def corranalysis(muts, cutoffnum, impsites, imp, filename, title):
    ## CORRELATION ANALYSIS - generates output png file correlation plot from information about sites, mutations, and the number of significant hits
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
    plt.savefig(f'{filename}.png', dpi=600);
    
    return corrs

def frame_cleaner(frame, sitenumlists, species, plist, betalist, anti):
    ## cleans a dataframe with data from several species analyses
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
    frame.to_csv(f'single_{anti}.csv', index = False)
    return frame

## Rifampicin Analyses
## all species analysis
# read in and clean data for analysis
muts, phens, rpobaln, rpobrefs = aln_cleaner("afas/rpob.afa", "combined_genomes.csv", "rifampicin")
# further clean and run logistic regression on data
allp, allmuts, allmutssp, allsites, allsitenum, allbetas = lg_regressor(rpobaln, muts, phens, rpobrefs, species = None, sp = 1)
# visualize the results
allrpob, pcut, allimpsites = framer(allp, allsites, allsitenum, rpobaln, rpobrefs, allbetas, 1, muts, "all", extimp = [])
qq_plot(allp, 'single_rif_all_qq', 'Multi-Species')
allcorrs = corranalysis(allmutssp, pcut, allimpsites, allrpob.Site, 'single_rif_all_corr', 'Multi-Species')

## tb analysis
# further clean and run logistic regression on data
tbp, tbmuts, tbmutssp, tbsites, tbsitenum, tbbetas = lg_regressor(rpobaln, muts, phens, rpobrefs, species = "M. tuberculosis")
# visualize results
tbrpob, tbpcut, tbimpsites = framer(tbp, tbsites, tbsitenum, rpobaln, rpobrefs, tbbetas, 1, tbmuts, "tb", extimp = list(allrpob.SiteNum))
newrpob = pd.merge(allrpob, tbrpob, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(tbp, 'single_rif_tb_qq', 'M. tuberculosis')
tbcorrs = corranalysis(tbmutssp, tbpcut, tbimpsites, tbrpob.Site, 'single_rif_tb_corr', 'M. tuberculosis')

## s. aureus analysis
# further clean and run logistic regression on data
stp, stmuts, stmutssp, stsites, stsitenum, stbetas = lg_regressor(rpobaln, muts, phens, rpobrefs, species = "S. aureus", refsp = 1)
# visualize results
strpob, stpcut, stimpsites = framer(stp, stsites, stsitenum, rpobaln, rpobrefs, stbetas, 3, stmuts, "staph", extimp = list(newrpob.SiteNum))
newrpob2 = pd.merge(newrpob, strpob, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(stp, 'single_rif_st_qq', 'S. aureus')
if stpcut > 1: 
    stcorrs = corranalysis(stmutssp, stpcut, stimpsites, strpob.Site, 'single_rif_st_corr', 'S. aureus')

## merge and clean tables
numall = []
numtb = []
numstaph = []
for i in newrpob2.SiteNum: 
    numall.append(sum(muts[:,i]))
    numtb.append(sum(muts[muts[:,-2] == 1][:,i]))
    numstaph.append(sum(muts[(muts[:,-2] == 0) & (muts[:,-3] == 0) & (muts[:,-1] == 0)][:,i]))
frame = pd.DataFrame(list(zip(newrpob2.Site, numall, numtb, numstaph)), 
                                columns = ['Site', 'num-all', 'num-tb', 'num-staph']) 
newrpob3 = pd.merge(newrpob2, frame, how = 'outer', on = ['Site'])
newrpob3 = newrpob3.astype({'num-all': 'int32', 'num-tb': 'int32', 'num-staph': 'int32'})
sitenumlists = [list(tbsitenum), list(stsitenum), list(allsitenums)]
species = ['tb', 'staph', 'all']
plist = [tbp, stp, allp]
betalist = [tbbetas, stbetas, allbetas]
anti = "rif"
newrpob4 = frame_cleaner(newrpob3, sitenumlists, species, plist, betalist, anti)


## Ciprofloxacin Analyses
## all species analysis
# read in and clean data for analysis
gyrmuts, gyrphens, gyraln, gyrarefs = aln_cleaner("afas/gyra.afa", "combined_genomes.csv", "ciprofloxacin")
# further clean and run logistic regression on data
allp, allmuts, allmutssp, allsites, allsitenums, allbetas = lg_regressor(gyraln, gyrmuts, gyrphens, gyrarefs, species = None, sp = 0, refsp = 0)
# visualize the results
allgyra, pcut, allimpsites = framer(allp, allsites, allsitenums, gyraln, gyrarefs, allbetas, 0, allmuts, "all", extimp = [])
qq_plot(allp, 'single_cip_all_qq', 'Multi-Species')
allcorrs = corranalysis(allmutssp, pcut, allimpsites, allgyra.Site, 'single_cip_all_corr', 'Multi-Species')

## e. coli analysis
# further clean and run logistic regression on data
tbp, tbmuts, tbmutssp, tbsites, tbsitenum, tbbetas = lg_regressor(gyraln, gyrmuts, gyrphens, gyrarefs, species = "E. coli", refsp = 0)
# visualize results
tbrpob, tbpcut, tbimpsites = framer(tbp, tbsites, tbsitenum, gyraln, gyrarefs, tbbetas, 0, tbmuts, "ec")
newrpob = pd.merge(allgyra, tbrpob, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(tbp, 'single_cip_ec_qq', 'E. coli')
tbcorrs = corranalysis(tbmutssp, tbpcut, tbimpsites, tbrpob.Site, 'single_cip_ec_corr', 'E. coli')

## s. enterica analysis
# further clean and run logistic regression on data
stp, stmuts, stmutssp, stsites, stsitenum, stbetas = lg_regressor(gyraln, gyrmuts, gyrphens, gyrarefs, species = "S. enterica", refsp = 0)
# visualize results
strpob, stpcut, stimpsites = framer(stp, stsites, stsitenum, gyraln, gyrarefs, stbetas, 0, stmuts, "sal")
newrpob2 = pd.merge(newrpob, strpob, how = 'outer', on = ['Site', 'SiteNum', 'RegAA'])
qq_plot(stp, 'single_cip_sal_qq', 'S. enterica')
if stpcut > 1: 
    stcorrs = corranalysis(stmutssp, stpcut, stimpsites, strpob.Site, 'single_cip_sal_corr', 'S. enterica')

## merge and clean tables
# get counts of mutants at each site
numall = []
numec = []
numsal = []
for i in newrpob2.SiteNum: 
    numall.append(sum(allmuts[:,i]))
    numec.append(sum(allmuts[allmuts[:,-3] == 1][:,i]))
    numsal.append(sum(allmuts[allmuts[:,-1] == 1][:,i]))
# add counts to dataframes, change type
frame = pd.DataFrame(list(zip(newrpob2.SiteNum, numall, numec, numsal)), 
                                columns = ['SiteNum', 'num-all', 'num-ec', 'num-sal']) 
newrpob3 = pd.merge(newrpob2, frame, how = 'outer', on = ['SiteNum'])
newrpob3 = newrpob3.astype({'num-all': 'int32', 'num-ec': 'int32', 'num-sal': 'int32'})
# fill in remaining blanks in merged table, where appropriate
sitenumlists = [list(tbsitenum), list(stsitenum), list(allsitenums)]
species = ['ec', 'sal', 'all']
plist = [tbp, stp, allp]
betalist = [tbbetas, stbetas, allbetas]
anti = "cip"
newrpob4 = frame_cleaner(newrpob3, sitenumlists, species, plist, betalist, anti)
