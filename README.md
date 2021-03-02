# STEAP 
**S**ingle cell **T**ype **E**nrichment **A**nalysis for **P**henotypes (**STEAP**) is an extension to [CELLECT](https://github.com/perslab/CELLECT). 
It contains additional post-processing scripts and uses an updated MAGMA version. 
CELLECT uses MAGMA v1.07b which is known to inflate false positives.
STEAP updates to MAGMA v1.08 which resolved these issues.

[de Leeuw, C., Sey, N. Y. A., Posthuma, D., & Won, H. (2020). A response to Yurko et al: H-MAGMA, inheriting a shaky statistical foundation, yields excess false positives. Cold Spring Harbor Laboratory. https://doi.org/10.1101/2020.09.25.310722](https://www.biorxiv.org/content/10.1101/2020.09.25.310722v1)

STEAP runs multiple post-processing steps on top of CELLECT output files:
  - GSEA
  - Cell-Type Correlation
  - ES Gene Correlation

![pipeline](https://github.com/erwinerdem/STEAP/blob/master/pipeline.png)


## Installation
**Step 1: Install (mini)conda**  
The pipeline works in conda enviroments, so you will have to install anaconda or miniconda first if not already on the system.
```
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh"
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh
```

Create a new conda enviroment
```
conda create -n steap
```



**Step 2: Clone CELLECT repository**  
Clone the repository: 
```
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
The `--recurse-submodules` is needed to clone the [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) 'ldsc' ([pascaltimshel/ldsc](https://github.com/pascaltimshel/ldsc)), which is a modfied version of the original ldsc repository.
(Cloning the repo might take few minutes as the CELLECT data files (> 1-3 GB) will be downloaded. To skip downloading the data files, use `GIT_LFS_SKIP_SMUDGE=1 git clone --recurse-submodules https://github.com/perslab/CELLECT.git` instead.)

...

## Getting Started
...
### Munging GWAS
...
### Converting scRNA Data to ES Matrix
Example notebooks converting the raw scRNA to ES matrices using [CELLEX](https://github.com/perslab/CELLEX) can be found in the [cellex-notebooks](https://github.com/erwinerdem/cellex-notebooks).

