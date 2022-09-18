#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

reference_genome=/iblm/netapp/data1/jezhou/cellranger/rn6
fastqs_dir=/iblm/netapp/data1/jezhou/Telese_Rat_Amygdala/snATAC/fastq
outdir=/iblm/netapp/data1/jezhou/Telese_Rat_Amygdala/snATAC/realigned_outputs

rat=Rat_Amygdala__cocaine_low
sample_names=Rat_Amygdala__cocaine_low

echo "ID=$rat"
echo "sample=$sample_names"

cd $outdir

cellranger-atac count --id=$rat \
--reference=$reference_genome \
--fastqs=$fastqs_dir \
--sample=$sample_names \
--localcores=8 \
--localmem=64 

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "realigning snATAC-seq FASTQs for $rat finished successfully" &>/dev/null
	else
		slack "realigning snATAC-seq FASTQs for $rat exited with error code $exit_code"
	fi
fi
exit "$exit_code"