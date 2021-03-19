# STEAP 
**S**ingle cell **T**ype **E**nrichment **A**nalysis for **P**henotypes (**STEAP**) uses scRNA-seq data and GWAS summary statistics to determine which cell-types are enriched in the GWAS phenotype. It is an extension to [CELLECT](https://github.com/perslab/CELLECT) and uses [S-LDSC](https://github.com/bulik/ldsc) ([Finucane et al., 2015](https://www.nature.com/articles/ng.3404)), [MAGMA](https://ctg.cncr.nl/software/magma) [(de Leeuw et al., 2015)](https://doi.org/10.1371/journal.pcbi.1004219) and [H-MAGMA](https://github.com/thewonlab/H-MAGMA) [(Sey et al., 2020)](https://doi.org/10.1038/s41593-020-0603-0) for enrichment analysis.


STEAP runs multiple post-processing steps on top of CELLECT output files:
  - Gene Set Enrichment Analysis (GSEA)
  - Cell-Type Correlation
  - Expression Specificity (ES) Gene Correlation

![pipeline](https://github.com/erwinerdem/STEAP/blob/master/pipeline.png)


## Installation
**Step 1: Install (mini)conda**  
The pipeline works in conda enviroments, so you will have to install anaconda or miniconda first if not already on the system.
You can do this by running:
```
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh"
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh
```

**Step 2: Clone this repo and run the install script** 

First clone this repo, and then run the install.sh script by running:
```
git clone https://github.com/erwinerdem/STEAP.git
bash STEAP/install.sh
conda activate steap
```
This will create  the necessary conda environments and clones the CELLECT repo and merges it with STEAP.

## Getting Started
Any GWAS summary statistic of your interest can be used in this tutorial.
Here, we will download and use the PGC depression GWAS as an example. Download can take a few minutes.
```
wget -O gwas/PGC_UKB_depression.txt https://datashare.is.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt --no-check-certificate
```

#### Munging GWAS
The summary statistic needs to be compatible with the pipeline. This is done by first munging the GWAS using the [`mtag_munge.py`](https://github.com/pascaltimshel/ldsc/blob/d869cfd1e9fe1abc03b65c00b8a672bd530d0617/mtag_munge.py) script from ldsc. This must be done in a special conda environment.
```
conda activate munge_ldsc
```
Then we can munge the PGC depression GWAS using:
```
python2 ldsc/mtag_munge.py \
--sumstats gwas/PGC_UKB_depression.txt \
--merge-alleles data/ldsc/w_hm3.snplist \
--a1 A1 \
--a2 A2 \
--snp MarkerName \
--p P \
--N-cas 170756 \
--N-con 329443 \
--signed-sumstats LogOR,0 \
--frq Freq \
--out gwas/PGC_UKB_depression
```
Each GWAS has different column names so change the parameters accordingly. More examples can be found in [timshel-2020/src/prep-gwas_munge/README-cmds_munge_sumstats.txt](https://github.com/perslab/timshel-2020/blob/master/src/prep-gwas_munge/README-cmds_munge_sumstats.txt). The [`mtag_munge.py`](https://github.com/pascaltimshel/ldsc/blob/d869cfd1e9fe1abc03b65c00b8a672bd530d0617/mtag_munge.py) script is also well documented.
The output `.sumstats.gz` file will be used as input for the pipeline.
Deactivating the munge_ldsc environment can be done simply using 
`conda deactivate`.

#### Converting scRNA Data to ES Matrix
Converting the scRNA to ES matrix requires CELLEX. This can be installed as mentioned in the [CELLEX repository](https://github.com/perslab/CELLEX#setup). Example notebooks converting the raw scRNA to ES matrices using CELLEX can be found in the [cellex-notebooks](https://github.com/erwinerdem/cellex-notebooks).

#### Running the pipeline
To run the pipeline a `config.yml` must first be set up. For the PGC depression GWAS we will use the [config.yml](https://github.com/erwinerdem/STEAP/tree/master/config/config.yml) file. This file can be edited to instead include your own GWAS summary statistics or scRNA-seq data.
Under `SPECIFICITY_INPUT` you can add the `id` and `path` of your ES matrix file. The `id` is simply the identifier for the dataset and the `path` the path to the file. 
Similarly, under `GWAS_SUMSTATS` you can add new munged GWAS summary statistics.
If you want to include a different annotation file for H-MAGMA, you can edit this under `HMAGMA_ANNOT`. The results will be saved in the `BASE_OUTPUT_DIR`.
For normal use it is recommended to only edit these parameters and leave the rest as is, unless you know what you are doing.

To run the enrichment analysis run:
```
snakemake --use-conda -j -s cellect-magma.snakefile --configfile config/config.yml
snakemake --use-conda -j -s cellect-h-magma.snakefile --configfile config/config.yml
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config/config.yml
```

#### Running the post-processing notebook
To use the post-processing scripts just use the [notebook](https://github.com/erwinerdem/STEAP/blob/master/notebooks/depression_example.ipynb) after running the pipeline.
