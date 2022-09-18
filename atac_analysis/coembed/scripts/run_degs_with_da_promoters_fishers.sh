#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-13
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $HOME/rat_snatacseq_analysis/coembed/celltypes.txt`
degs=$HOME/rat_snrnaseq_pipeline/dge_labeled/out/dge-with-batch-percent-sample-covars/MAST/high_vs_low/${celltype}_high_vs_low_with_fdr.csv
dapeaks=$HOME/rat_snatacseq_analysis/dapeak_analysis_rn6/out/full_cocaine_dataset_with_naive_rn6_with_sample_covar_prefiltered/da_peaks/high_vs_low/negbinom/${celltype}_dapeaks_with_fdr.csv

outdir=$HOME/rat_snatacseq_analysis/coembed/fishers_degs_with_da_promoters_with_sample_covar
mkdir -p $outdir
outfh=$outdir/${celltype}_fishers.csv

echo "running fishers test measuring enrichment of DEGs with DA promoters in $celltype"

Rscript degs_with_da_promoters_fishers.R --celltype $celltype \
--degs $degs --dapeaks $dapeaks --outfh $outfh

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "fishers exact test finished successfully for $celltype (enrichment of DEGs with DA promoters)" &>/dev/null
    else
        slack "fishers exact test exited error code $exit_code for $celltype (enrichment of DEGs with DA promoters)"
    fi
fi
exit "$exit_code"
