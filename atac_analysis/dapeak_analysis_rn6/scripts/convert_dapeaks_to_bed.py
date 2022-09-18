import argparse
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

parser = argparse.ArgumentParser()
parser.add_argument("--dapeaks", action = "store",
	help = "da peaks csv file to be converted to bed")
parser.add_argument("--celltype", action = "store",
	help = "cell type")
parser.add_argument("--bed", action = "store",
	help = "where to save output bed")
parser.add_argument("--significant", action = "store_true",
	help = "only convert significant peaks")

args = parser.parse_args()

df = pd.read_csv(args.dapeaks, index_col = 0)

if args.significant:
	res = fdrcorrection(df['p_val'])
	df = df[res[0]]

pk_ranges = [['chr' + str(x[0]), x[1], x[2]] for x in [x.split('-') for x in df.index]]

pk_names = [("%s_peak_%d" % (args.celltype, x)) for x in list(range(1,df.shape[0]+1))]

bed = pd.DataFrame(pk_ranges)

bed['name'] = pk_names

bed.to_csv(args.bed, header = False, index = False, sep = "\t")