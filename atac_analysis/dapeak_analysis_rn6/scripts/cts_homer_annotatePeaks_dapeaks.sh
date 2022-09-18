#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-12
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snatacseq_analysis/dapeak_analysis_rn6/out/full_cocaine_dataset_with_naive_rn6_with_sample_covar-X/da_peaks/high_vs_low/negbinom
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $HOME/rat_snatacseq_analysis/dapeak_analysis_rn6/data/rds/full_cocaine_dataset_with_naive_rn6/celltypes_list.txt`
infile=$wd/${celltype}_dapeaks_significant.csv
gtf=/iblm/netapp/data1/jezhou/cellranger/rn6-2014-build/Rattus_norvegicus.Rnor_6.0.98.gtf.filtered
output=$wd/${celltype}_significant_homer_annotatePeaks_output.txt

echo $infile

cut -f1 -d, $infile | tr '-' '\t' > $wd/${celltype}_significant.bed && \
perl $HOME/homer/bin/annotatePeaks.pl $wd/${celltype}_significant.bed rn6 \
-gtf $gtf > $output

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "running HOMER annotatePeaks.pl for DA peaks finished successfully for $celltype" &>/dev/null
	else
		slack "running HOMER annotatePeaks.pl for DA peaks exited with error code $exit_code for $celltype"
	fi
fi
exit "$exit_code"
