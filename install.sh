#!/bin/sh
git lfs install 
git lfs clone --recurse-submodules https://github.com/perslab/CELLECT.git 
rsync -a --remove-source-files STEAP/* CELLECT
rm -rf STEAP 
mv CELLECT/ STEAP/
cd STEAP
conda env create -f environment_steap.yml # create conda env needed for running STEAP
conda env create -f ldsc/environment_munge_ldsc.yml # create the conda env for GWAS munging
