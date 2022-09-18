#!/bin/bash
#$ -q iblm.q
#$ -V
#$ -cwd

Rscript ./dge_ambiguous_clusters.R

# message the user on slack if possible
exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
		slack "pipeline finished successfully" &>/dev/null
	else
		slack "pipeline exited with error code $exit_code"
	fi
fi
exit "$exit_code"
