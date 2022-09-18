import argparse
import pandas as pd
from pyfaidx import Fasta
from statsmodels.stats.multitest import fdrcorrection
import sys

def get_sequence(region):
    print(region)
    key, start, end = region.split("-")
    tmpseq = chrs[key][int(start):int(end)].seq
    return tmpseq.upper()

parser = argparse.ArgumentParser()
parser.add_argument("--coords", action = "store", help = 'BED format file with coordinates of sequences to get')
parser.add_argument("--fasta", action = 'store', help = 'path to FASTA for genome')
# parser.add_argument("--bp", action = "store", type = int, 
# 	help = "window size for calculating cutoff space from ends of chromosomes")
parser.add_argument("--out", action = "store", help = "output file")
parser.add_argument("--bg", action = "store_true", 
    help = "if true, process differently (input files are formatted differently b/w fg/bg")
parser.add_argument("--alpha", action = "store", default = 0.1,
    help = "cutoff for selecting significant fg peaks")
parser.add_argument("--pos", action = "store_true", 
    help = "select positive fg peaks (logFC>0)")
parser.add_argument("--neg", action = "store_true",
    help = "select negative fg peaks (logFC<0")

args = parser.parse_args()

if args.pos and args.neg:
    print("cannot set both pos and neg flags as True")
    sys.exit()

if args.bg and args.pos:
    print("cannot set both bg and pos flag as True")
    sys.exit()

if args.bg and args.neg:
    print("cannot set both bg and neg flag as True")
    sys.exit()

print("loading FASTA")
chrs = Fasta(args.fasta)

if args.bg:
    """
    pick bg peaks 
    get seqs based on string 
    e.g. chr1-1-1000
    """
    print("getting bg seqs")
    regions = open(args.coords).read().splitlines()
    seqs = map(get_sequence, regions)
    with open(args.out, 'w') as fh:
        for peak,seq in zip(regions,seqs):
            fh.write(">" + peak + "\n" + seq + "\n")
else:
    """
    pick fg peaks
    get seqs from strings in differential analysis outputs
    """
    print("getting fg peaks")
    peaks = pd.read_csv(args.coords, index_col = 0)
    peaks.index.name = "regions"
    peaks.reset_index(inplace = True)
    # calculate q-val
    q_bool, q_val = fdrcorrection(peaks['p_val'], alpha=args.alpha, method='indep', is_sorted=False)
    peaks['q_val'] = q_val
    # select significant peaks (default: FDR<10%)
    significant = peaks[peaks['q_val']<args.alpha]
    if args.pos:
        print("getting positive fg peaks")
        significant = significant[significant['avg_log2FC']>0]
    elif args.neg:
        print("getting negative fg peaks")
        significant = significant[significant['avg_log2FC']<0]
    else:
        print("getting fg peaks regardless of directionality of foldchange")
    # peaks['sequence'] = peaks['regions'].apply(lambda x: get_sequence(x))
    significant['sequence'] = significant['regions'].apply(lambda x: get_sequence(x))
    with open(args.out, 'w') as fh:
        print("writing sequences to fasta file %s" % args.out)
        for peak,seq in zip(significant['regions'].to_list(), significant['sequence'].to_list()):
            fh.write(">" + peak + "\n" + seq + "\n")
    # peaks['sequence'].to_csv(args.out, index = False, header = False)