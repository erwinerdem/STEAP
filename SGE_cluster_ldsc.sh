#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V

snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config/config.yml
