import argparse
from os import path

parser = argparse.ArgumentParser()
parser.add_argument("--bg", action = "store",
	help = "string to write for bg in ldcts")
parser.add_argument("--bg_dir", action = "store",
	help = "directory where bg ldscore files are located")
parser.add_argument("--celltypes", action = "store", nargs = "+",
	help = "list of cts to include in ldcts")
parser.add_argument("--cts_dir", action = "store",
	help = "directory where cts ldscore files are located")
parser.add_argument("--out", action = "store",
	help = "where to write ldcts file")

args = parser.parse_args()

bg = path.join(args.bg_dir, "%s." % args.bg)

with open(args.out, 'w') as fh:
    for celltype in args.celltypes:
        fg = path.join(args.cts_dir, "%s." % celltype)
        line = "%s\t%s,%s\n" % (celltype, fg, bg)
        fh.write(line)
