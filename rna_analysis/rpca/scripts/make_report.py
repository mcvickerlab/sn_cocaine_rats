import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--samples', action = "store", nargs='+',
	help='list of samples that were integrated')
parser.add_argument('--subset', action='store_true',
	help='indicate whether samples were subset to remove low quality cells')
parser.add_argument("--rpca", action = "store_true",
	help = "indicate whether RPCA was used for integrating data")
parser.add_argument("--test", action = "store",
	help = "name of DGE test")
parser.add_argument("--prefilt", action = "store_true",
	help = "indicate whether genes were prefiltered before DGE analysis")
parser.add_argument("--covars", action = "store", nargs = "+",
	help = "list of covars used for dge test")
parser.add_argument("--outfh", action = "store",
	help = "where to write report to")

args = parser.parse_args()

with open(args.outfh, "w") as fh:
	samples_str = ','.join(args.samples)
	fh.write("samples: %s" % samples_str)
	fh.write("subset: %s" % args.subset)
	fh.write("rpca: %s" % args.rpca)
	fh.write("test: %s" % args.test)
	fh.write("prefilt: %s" % args.prefilt)
	covars_str = ",".join(args.covars)
	fh.write("covars: %s" % covars_str)
fh.close()


