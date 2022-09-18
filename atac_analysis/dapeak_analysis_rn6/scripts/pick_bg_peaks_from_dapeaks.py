import argparse
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

parser = argparse.ArgumentParser()
parser.add_argument("--coords", action = "store", 
	help = 'table with dapeak results')
parser.add_argument("--out", action = "store", help = "output file")
parser.add_argument("--alpha", action = "store", default = 0.1, 
    help = "cutoff for picking backgrounds")

args = parser.parse_args()

peaks = pd.read_csv(args.coords, index_col = 0)
peaks.index.name = "regions"
peaks.reset_index(inplace = True)
q_bool, q_val = fdrcorrection(peaks['p_val'], alpha=0.05, method='indep', is_sorted=False)
peaks['q_val'] = q_val

bg_peaks = peaks[peaks['q_val']>args.alpha]['regions'].tolist()

with open(args.out, 'w') as fh:
	fh.write('\n'.join(bg_peaks))