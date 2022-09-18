import pandas as pd
import glob
import os
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--sample", action = "store",
	help = "name of sample")
parser.add_argument("--results", action = "store", nargs = "+",
	help = "list of files with results to include in plot")
parser.add_argument("--outfh", action = "store",
	help = "where to save heatmap")
parser.add_argument("--ukb", action = "store_true",
    help = "if true, only plot results with UKB")
parser.add_argument("--log", action = "store_true",
    help = "if true, plot -log10(q) in heatmap")

args = parser.parse_args()
print(args.results)

studies = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in args.results]

dfs_list = []
for study,file in zip(studies, args.results):
    df = pd.read_csv(file, sep = '\t')
    df['Study'] = study
    fdrbh = fdrcorrection(df['Coefficient_P_value'], alpha = 0.1)
    df['fdr_correction'] = fdrbh[1]
    if args.ukb:
        if 'WithoutUKB' not in study:
            dfs_list.append(df)
        else: 
            continue
    else:
        dfs_list.append(df)

summary = pd.concat(dfs_list)
map_cts = {"Mature_Oligodendrocytes":"Oligodendrocytes",
"Precursor_Oligodendrocytes": "OPC"}
summary  = summary.replace({"Name":map_cts})

figwidth = 10
if args.sample == "pfc":
    figwidth = 14

if args.log:
    summary['-log10q'] = -np.log10(summary['fdr_correction'])
    plot_df = summary.pivot("Study", "Name", "-log10q")
    
    fig, ax1 = plt.subplots(1, 1, figsize = (figwidth,6), dpi=300)
    sns.heatmap(plot_df, cmap = "Blues", annot = True,
               cbar_kws={'label': '-log10(q)'})
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    ax1.tick_params(axis = 'x', labelrotation = 30)
    ax1.set_title(args.sample)
    plt.tight_layout()
    plt.savefig(args.outfh)
else:
    plot_df = summary.pivot("Study", "Name", "fdr_correction")

    fig, ax1 = plt.subplots(1, 1, figsize = (figwidth,6), dpi=300)
    sns.heatmap(plot_df, cmap = "Blues_r", annot = True,
               cbar_kws={'label': 'FDR-corrected p-value'})
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    ax1.tick_params(axis = 'x', labelrotation = 30)
    ax1.set_title(args.sample)
    plt.tight_layout()
    plt.savefig(args.outfh)
