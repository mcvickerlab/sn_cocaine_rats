#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd
#$ -t 1-11
#$ -tc 3
#$ -o /iblm/netapp/data1/jezhou/out
#$ -e /iblm/netapp/data1/jezhou/err

samples_table=$HOME/realign_snatac_fastqs/sample_names.txt
reference_genome=/iblm/netapp/data1/jezhou/cellranger/rn6
fastqs_dir=/iblm/netapp/data1/jezhou/Telese_Rat_Amygdala/snATAC/fastq
outdir=/iblm/netapp/data1/jezhou/Telese_Rat_Amygdala/snATAC/realigned_outputs

rat=`awk -v line=$SGE_TASK_ID 'NR==line' $samples_table | cut -f1`
sample_names=`awk -v line=$SGE_TASK_ID 'NR==line' $samples_table | cut -f2`

echo "ID=$rat"
echo "sample=$sample_names"

cd $outdir

# echo "cellranger-atac count --id=$rat \
# --reference=$reference_genome \
# --fastqs=$fastqs_dir \
# --sample=$sample_names \
# --localcores=8 \
# --localmem=64 >${rat}_log.out 2>${rat}_log.err"

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