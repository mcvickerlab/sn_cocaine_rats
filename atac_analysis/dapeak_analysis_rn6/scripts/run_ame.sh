#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-12
#$ -tc 4
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snatacseq_analysis/dapeak_analysis_rn6
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/data/rds/full_cocaine_dataset_with_naive_rn6/celltypes_list.txt`
outdir=${celltype}_ame_out_with_da_bg_vertebrates_no-thres
motifs=/iblm/netapp/data1/external/JASPAR2022/JASPAR2022_CORE_vertebrates_redundant_pfms_meme.txt

echo "running AME on $celltype"

cd $wd/out/full_cocaine_dataset_with_naive_rn6_no_sample_covar/da_peaks/high_vs_low/negbinom
ame --evalue-report-threshold 1500 --o $outdir \
--control ${celltype}-bg_seqs_from_dapeaks.fa ${celltype}-fg_seqs.fa $motifs

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "running AME finished successfully for $celltype" &>/dev/null
	else
		slack "running AME exited with error code $exit_code for $celltype"
	fi
fi
exit "$exit_code"
