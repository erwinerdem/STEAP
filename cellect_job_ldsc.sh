#!/bin/bash
# This is a sample submission script . Lines starting with # are
# comments . The first line ( with #!) should be in every script .
# Let 's set some variables for SGE. Lines starting with #$ are
# interpreted by SGE as if they were options to the qsub command
# (dn 't remove the # from the lines starting with #$).
#$ -S /bin/bash
#$ -cwd
#$ -M e.erdem@erasmusmc.nl
#$ -m be
#$ -j y
#$ -V

snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config/config.yml
