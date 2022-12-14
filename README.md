# MTLRank
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction
MTLRank is a multi-task learning based framework for inferring regulatory interactions from single cell data. MTLRank ranks TFs for each gene by learning models that predict RNA velocity values of target genes. Models are learned in a multi-task, soft-parameter-sharing based manner to improve the performance. MTLRank uses TF expression matrix and TF activity matrix as inputs and ranks the TFs with deep SHAP based on the trained models. This repository shows the Jupyter notebooks, and other script files used to test different models used in the study.
![workflow](/img/model_github_page.png)

## How to use
### Clone the Repository
```bash
git clone https://github.com/alexQiSong/MTLRank.git
```
### Environment Set up
1. Install singularity in your OS. We use singularity image to run the search of potential target genes by ChIPseeker. [Install singularity here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).
2. All required dependencies are specified in the file `MTLRANK_env.yml` in this repository. Download this file and use `conda` to recreate this environment by
```shell
conda env create -f MTLRank_env.yml
```
### Run pipeline
1. Activate the MTLRank environment by
```shell
conda activate MTLRank
```

2. Download sample data. We provide a script for downloading sample data sets of spleen and liver tissue from HuBMAP data portal. Due to data sharing policy, we are unable to directly share the data from HuBMAP consortium. Some of the data sets we used in the paper are not yet published. So the sample data sets might produce results different than what have been presented in the paper. To download the sample data sets, simply run `download.py` in this repository by (make sure MTLRank environment is active)
```shell
python download.py
```
This will automatically download most of the required data.

3. Download data from Cistrome database. Due to data sharing policy, we are not allowed to directly share ChIP-seq data from CistromeDB. Please go to http://cistrome.org/db/#/bdown and select Human_Factor to download all bed files and a QC file for ChIP-seq data. All bed files are downloaded as a gz file named "human_factor.tar.gz". Download and decompress this file and move all bed files to data/chipseq_bed/. QC file is named as "human_factor_full_QC.txt". Download this file and move it to data/chipseq_qc/

4. Run preprocessing steps. Simply run `preprocess.py` to perform all preprocessing steps including generating RPKM values and TF activity score computation. You may use `n_jobs` to perform parallel computation. This may take some time (~ 1hr when running with `--n_jobs 30`)
```shell
python preprocess.py --n_jobs 30
```

5. Run MTLRank pipeline using the jupyter notebookes `MTLRank_step1_train.ipynb`, `MTLRank_step2_rankTF.ipynb`, and `MTLRank_step3_GRN.ipynb` by following the order specified by the names. We provide more detailed instructions in each notebook, including the the suggested settings of hyperparameters.

### Contact
Contact us if you have any questions:  
Qi (Alex) Song: qisong@andrew.cmu.edu; sqsq3178@gmail.com  
Ziv Bar-Joseph: zivbj@andrew.cmu.edu  

### Cite
Comming soon!

### Copyright
©2022 Qi Song, Ziv Bar-Joseph. [Systems Biology Group at Carnegie Mellon University](http://www.sb.cs.cmu.edu/). All rights reserved.
