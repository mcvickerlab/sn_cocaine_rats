#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-10
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snatacseq_analysis/coembed
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/fimo_celltypes.txt`
input=${wd}/fimo/${celltype}_peak_seqs/da_peak_seqs.fa
outfh=${wd}/fimo/${celltype}_peak_seqs/da_peak_seqs_markov_bg.txt

fasta-get-markov -m 1 $input $outfh

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "fasta-get-markov finished successfully for $celltype (DA peaks linked to a DEG)" &>/dev/null
	else
		slack "fasta-get-markov exited with error code $exit_code for $celltype (DA peaks linked to a DEG)"
	fi
fi
exit "$exit_code"