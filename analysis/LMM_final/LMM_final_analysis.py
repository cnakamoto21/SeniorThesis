## IMPORT RELEVANT PACKAGES / FUNCTIONS
from evcouplings.align import Alignment, map_matrix
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv

## generate png qq plot from p-value data and labels
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

## generate the dataframe that displays the results for a single analysis
def framer(pvals, ORs, sites, sitenum, group):
    cutoffnum = sum([1 if i < (0.05 / len(pvals)) else 0 for i in sorted(pvals)])
    impsites = np.argsort(pvals).tolist()
    impnames = [sites[i] for i in impsites[0:cutoffnum]]
    impnums = [sitenum[i] for i in impsites[0:cutoffnum] if sitenum[i]]
    ps = sorted(pvals)[0:cutoffnum]
    ORs = list(np.array(ORs)[np.argsort(pvals)])
    frame = pd.DataFrame(list(zip(impnames, impnums, ORs, ps)), 
                                columns = ['Site', 'SiteNum', f'OR-{group}', f'p-{group}']) 
    return frame, cutoffnum, impsites

## generate png correlation plot with mutation data, labels
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

## data frame cleaning function: adds in missing data where possible, formats
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
    frame.to_csv(f'wgs_{anti}.csv', index = False)
    return frame

## data frame cleaning function for false positive analysis: adds missing data, formats
def frame_cleaner_fp(frame, sitenumlists, species, plist, betalist, anti):
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
    frame.to_csv(f'wgs_fp_{anti}.csv', index = False)
    return frame

## call functions
frames = []
altspclist = [['all', 'M. tb', 'S. aureus'], ['all', 'E. coli', 'S. enterica']]
superlist = [['all', 'M.tb', 'S.aureus'], ['all', 'E.coli', 'S.enterica']]
# iterate over antibiotics, 
for idx, anti in enumerate(['rifampicin', 'ciprofloxacin']): 
    plist, betalist, sitenums, counts = [], [], [], []
    levellist = superlist[idx]
    for i, species in enumerate(levellist): 
        muts = np.loadtxt(f'gemma_inputs_final/{species}_{anti}.loci', 'str', delimiter = ',')
        mutpd = pd.DataFrame(muts.T[3:,:], columns = list(muts.T[0,:]))
        sitepd = pd.read_csv(f'gemma_int_outputs/{anti[0:3]}_genedrop_{species}_siteinfo.csv')
        out = pd.read_csv(f'gemma_final_outputs/{species}_{anti}_gemma_out.assoc.txt', delimiter = '\t')
        sitenumtemp = [int(sitepd.loc[sitepd.Site == a].SiteNum) if len(sitepd.loc[sitepd.Site == a]) == 1
               else None for a in out.rs]
        newframe, cutoffnum, impsites = framer(out.p_wald, np.exp(np.array(out.beta)), out.rs, sitenumtemp, 
                                               species)
        qq_plot(out.p_wald, species, f'wgs_{anti}_{species}_qq')
        if cutoffnum > 1: 
            corranalysis(mutpd, cutoffnum, impsites, newframe.Site, species, f'wgs_{anti}_{species}_corr')
        plist.append(out.p_wald)
        betalist.append(np.exp(np.array(out.beta)))
        sitenums.append(sitenumtemp)
        counts.append(sitepd)
        if i > 0: 
            frame = pd.merge(frame, newframe, how = 'outer', on = ['Site', 'SiteNum'])
        else: 
            frame = newframe
    for i, species in enumerate(altspclist[idx]): 
        if species == "M. tb": species = "M. tuberculosis"
        frame[f'num-{levellist[i]}'] = [int(counts[i].loc[counts[i].Site == r.Site][f'num {species}']) if 
                                        len(counts[i].loc[counts[i].Site == r.Site]) == 1 else None
                                        for index,r in frame.iterrows()]
    frame = frame.astype({f'num-{levellist[i]}': 'int32'})
    frame = frame_cleaner(frame, sitenums, levellist, plist, betalist, anti)
    frames.append(frame)
