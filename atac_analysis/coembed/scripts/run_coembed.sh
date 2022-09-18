#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

atac=$HOME/rat_snatacseq_analysis/get_combined_peak_mtx/data/rds/full_cocaine_dataset_with_naive_rn6/harmony_with_predicted_id.rds
rna=$HOME/rat_snrnaseq_pipeline/rpca/out/subset/integrated_relabeled_inhneuron_subtypes.rds
out=$HOME/rat_snatacseq_analysis/get_combined_peak_mtx/data/rds/full_cocaine_dataset_with_naive_rn6/coembedding.rds

Rscript coembed.R --atac $atac --rna $rna --out $out --activity

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "coembedding finished successfully" &>/dev/null
	else
		slack "coembedding exited with error code $exit_code"
	fi
fi
exit "$exit_code"