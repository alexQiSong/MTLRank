# MTLRank
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction
MTLRank is a multi-task learning based framework for inferring regulatory interactions from single cell data. MTLRank ranks TFs for each gene by learning models that predict RNA velocity values of target genes. Models are learned in a multi-task, soft-parameter-sharing based manner to improve the performance. MTLRank uses TF expression matrix and TF activity matrix as inputs and ranks the TFs with deep SHAP based on the trained models. This repository shows the Jupyter notebooks, and other script files used to test different models used in the study.

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
2.  
