#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-13
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snatacseq_analysis/coembed
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/data/rds/full_cocaine_dataset_with_naive_rn6/celltypes_list.txt`
rds=$wd/data/rds/full_cocaine_dataset_with_naive_rn6/coembedding.rds
outdir=$wd/peak_to_gene_links

mkdir -p $outdir

Rscript link_peaks_to_genes_cts.R --rds $rds --celltype $celltype --out $outdir/$celltype.rds

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "make peak bc matrix from combined finished successfully for $celltype" &>/dev/null
	else
		slack "make peak bc matrix from combined exited with error code $exit_code for $celltype"
	fi
fi
exit "$exit_code"