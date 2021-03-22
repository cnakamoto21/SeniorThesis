## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv

## SOME UNIVERSAL STUFF!
# make a dictionary encoding the order of the species
sptoid = {"562": 0, "1773": 1, "1733": 1, "28901": 2, "1280": 3}
sptoname = {"562": "E. coli", "1773": "M. tb", "1733": "M. tb", "28901": "S. enterica", "1280": "S. aureus"}
speciestosp = {"S. aureus": 3, "E. coli": 0, "M. tuberculosis": 1, "S. enterica": 2}

## make qq plot png file based on p-value data
def qq_plot(p_val_list, title, filename, xlim=[0,10]):
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

## make dataframe with results for a given level of control, antibiotic
def framer_alt(pvals, ORs, sites, sitenum, group):
    cutoffnum = sum([1 if i < (0.05 / len(pvals)) else 0 for i in sorted(pvals)])
    impsites = np.argsort(pvals).tolist()
    impnames = [sites[i] for i in impsites[0:cutoffnum]]
    impnums = [sitenum[i] for i in impsites[0:cutoffnum] if sitenum[i]]
    ps = sorted(pvals)[0:cutoffnum]
    ORs = list(np.array(ORs)[np.argsort(pvals)])
    frame = pd.DataFrame(list(zip(impnames, impnums, ORs, ps)), 
                                columns = ['Site', 'SiteNum', f'OR-{group}', f'p-{group}']) 
    return frame, cutoffnum, impsites

## make png correlation plot for given mutation, p-value cutoff data and labels
def corranalysis(muts, cutoffnum, imp, sitenames, title, filename):
    ## CORRELATION ANALYSIS
    # note: uses variable definitions from all species set
    importantsites = np.empty([cutoffnum, muts.shape[0]])
    for i in range(cutoffnum):
        importantsites[i,:] = muts.iloc[:,[imp[i]]].values.reshape(-1,)
    corrs = np.corrcoef(importantsites)

    import seaborn as sns
    import matplotlib.pyplot as plt

    # plot the heatmap
    fig = plt.figure(figsize=(11,10))
    ax = fig.gca
    sns.heatmap(corrs, cmap="Blues", xticklabels=sitenames, yticklabels=sitenames).set_title(f'{title}')
    plt.savefig(f'figures/{filename}.png', dpi=600);
    
    return corrs

## clean merged frame: format data, get insignificant ORs and p-values
def frame_cleaner_alt(frame, sitenumlists, species, plist, betalist, anti):
    for index, row in frame.iterrows():
        for i in range(len(species)):
            if np.isnan(row[f'p-{species[i]}']):
                if row['SiteNum'] in sitenumlists[i]: 
                    frame.loc[index, f'p-{species[i]}'] = "{:.1E}".format(plist[i][sitenumlists[i].index(row['SiteNum'])])
                    frame.loc[index, f'OR-{species[i]}'] = "{:.2f}".format(betalist[i][sitenumlists[i].index(row['SiteNum'])])
            elif isinstance(row[f'p-{species[i]}'], float): 
                frame.loc[index, f'p-{species[i]}'] = "{:.1E}*".format(frame.loc[index, f'p-{species[i]}'])
                frame.loc[index, f'OR-{species[i]}'] = "{:.2f}".format(frame.loc[index, f'OR-{species[i]}'])
    collist = ['Site', 'num all'] + [f'OR-{s}' for s in species] + [f'p-{s}' for s in species]
    frame = frame[collist]
    frame.to_csv(f'inter_{anti}.csv', index = False)
    return frame

## iteratively call functions to generate figures and table
frames = []
levellist = ['genedrop', 'classic', 'dropsite']
for anti in ['rif', 'cip']:
    # get binary encoded matrix
    muts = np.loadtxt(f'gemma_int_outputs/{anti}_inter_all.loci', 'str', delimiter = ',')
    mutpd = pd.DataFrame(muts.T[3:,:], columns = list(muts.T[0,:]))
    # get information about the sites of interest in the data
    sitepd = pd.read_csv(f'gemma_int_outputs/{anti}_all_siteinfo.csv', index_col=0)
    plist, betalist, sitenums = [], [], []
    for level in levellist:
        # get significant hits
        out = pd.read_csv(f'gemma_int_outputs/{anti}_{level}_all_gemma_out.assoc.txt', delimiter = '\t')
        # call functions for table, plots for the analysis
        newframe, cutoffnum, impsites = framer_alt(out.p_wald, np.exp(np.array(out.beta)), out.rs, 
                                                   sitepd.SiteNum, level)
        qq_plot(out.p_wald, level, f'inter_{anti}_{level}_qq')
        corranalysis(mutpd, cutoffnum, impsites, newframe.Site, level, f'inter_{anti}_{level}_corr')
        # store relevant information for big table
        plist.append(out.p_wald)
        betalist.append(np.exp(np.array(out.beta)))
        sitenums.append(list(sitepd.SiteNum))
        # merge table with others if appropriate
        if level != 'genedrop': 
            frame = pd.merge(frame, newframe, how = 'outer', on = ['Site', 'SiteNum'])
        else: 
            frame = newframe
    # get the counts, add to dataframe
    frame['num all'] = [sitepd.loc[sitepd.Site == r.Site]['num all'] for i,r in frame.iterrows()]
    frame = frame.astype({'num all': 'int32'})
    # call frame cleaner to format, incorporate insignificant p-values and ORs
    frame = frame_cleaner_alt(frame, sitenums, levellist, plist, betalist, anti)
    frames.append(frame)
