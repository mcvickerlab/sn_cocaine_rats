import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gex", action = 'store', help = 'sample name for GEX file')
parser.add_argument("--atac", action = "store", help = "sample name for ATAC file")
parser.add_argument("--fastqs", action = "store", help = "path to directory containing FASTQs")
parser.add_argument("--out", action = "store", help = "output file name")

args = parser.parse_args()
print(args)


with open(args.out, 'w') as fh:
    fh.write('fastqs,sample,library_type\n')
    fh.write(args.fastqs + ',' + args.gex + ',Gene Expression\n')
    fh.write(args.fastqs + ',' + args.atac + ',Chromatin Accessibility\n')

