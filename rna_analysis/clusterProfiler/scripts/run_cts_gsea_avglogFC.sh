#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd

outdir=$HOME/rat_snrnaseq_pipeline/clusterProfiler/dge-with-batch-percent-sample-covars/MAST/high_vs_low/gsea_avg_logFC/results

mkdir -p $outdir

for f in $HOME/rat_snrnaseq_pipeline/dge_labeled/out/dge-with-batch-percent-sample-covars/MAST/high_vs_low/*.csv; do
	celltype=$(basename $f _high_vs_low.csv)
	echo "running GSEA on $celltype"
	Rscript $PWD/gsea.R --markers $f --kegg --go --cutoff 1 --outdir $outdir --celltype $celltype 
done

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "running cell type-specific GWAS based on avg_logFC finished successfully" &>/dev/null
	else
		slack "running cell type-specific GWAS based on avg_logFC exited with error code $exit_code"
	fi
fi
exit "$exit_code"
