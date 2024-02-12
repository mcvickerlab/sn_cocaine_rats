import pandas as pd
import glob
import os
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample", action = "store",
	help = "name of sample")
parser.add_argument("--results", action = "store", nargs = "+",
	help = "list of files with results to include in plot")
parser.add_argument("--outfh", action = "store",
	help = "where to save heatmap")
parser.add_argument("--ukb", action = "store_true",
    help = "if true, only plot results with UKB")

args = parser.parse_args()


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
plot_df = summary.pivot("Study", "Name", "fdr_correction")

p = sns.clustermap(plot_df, row_cluster = False, cmap = sns.cm.rocket_r,
           cbar_kws={'label': 'FDR corrected p-value'})
p.fig.set_size_inches(12,6)
p.ax_heatmap.set_xlabel('')
p.ax_heatmap.set_ylabel('')
p.ax_heatmap.tick_params(axis = 'x', labelrotation = 30)
p.fig.suptitle(args.sample)
p.savefig(args.outfh)