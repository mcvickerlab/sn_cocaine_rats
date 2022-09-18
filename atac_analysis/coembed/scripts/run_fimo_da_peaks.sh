#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-10
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

wd=$HOME/rat_snatacseq_analysis/coembed
celltype=`awk -v line=$SGE_TASK_ID 'NR==line' $wd/fimo_celltypes.txt`
bfile=${wd}/fimo/${celltype}_peak_seqs/da_peak_seqs_markov_bg.txt
# motifs=${wd}/fimo/${celltype}_motifs_meme.txt
motifs="/iblm/netapp/data1/external/JASPAR2022/JASPAR2022_CORE_vertebrates_redundant_pfms_meme.txt"
sequence=${wd}/fimo/${celltype}_peak_seqs/da_peak_seqs.fa
outdir=${wd}/fimo/${celltype}_fimo_out_da_peaks_all_JASPAR_motifs

fimo --bfile $bfile --o $outdir $motifs $sequence

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "fimo finished successfully for $celltype (DA peaks linked to a DEG, all JASPAR motifs)" &>/dev/null
	else
		slack "fimo with error code $exit_code for $celltype (DA peaks linked to a DEG, all JASPAR motifs)"
	fi
fi
exit "$exit_code"