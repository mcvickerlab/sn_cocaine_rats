#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-13
#$ -tc 3
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snrnaseq_pipeline/dge_labeled
outdir=$wd/out/relabeled_subtypes3
rds=$HOME/rat_snrnaseq_pipeline/rpca/out/subset/integrated_relabeled_inhneuron_subtypes_with_batch.rds
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/celltypes_list.txt`
# mast_dir=/iblm/netapp/home/jezhou/rat_snrnaseq_pipeline/dge_labeled/out/relabeled_subtypes/high_vs_low
mkdir -p $outdir

echo "Running permutation test for DEGs in $celltype"

Rscript dge_celltypes_cts.R --celltype $celltype --outdir $outdir --rds $rds

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "mast_permutation_test.R finished successfully for $celltype" &>/dev/null
	else
		slack "mast_permutation_test.R exited with error code $exit_code for $celltype"
	fi
fi
exit "$exit_code"
