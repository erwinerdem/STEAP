#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V

snakemake --use-conda -j -s cellect-h-magma.snakefile --configfile config/config.yml
