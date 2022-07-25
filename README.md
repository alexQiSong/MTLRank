# MTLRank
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction
MTLRank is a multi-task learning based framework for inferring regulatory interactions from single cell data. MTLRank ranks TFs for each gene by learning models that predict RNA velocity values of target genes. Models are learned in a multi-task, soft-parameter-sharing based manner to improve the performance. MTLRank uses TF expression matrix and TF activity matrix as inputs and ranks the TFs with deep SHAP based on the trained models. This repository shows the Jupyter notebooks, and other script files used to test different models used in the study.

## How to use
### Clone the Repository
```bash
git clone https://github.com/alexQiSong/MTLRank.git
```
### R Environment Set up
1. Please install R package "ChIPSeeker". More information can be found [here](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#:~:text=ChIPseeker%20is%20an%20R%20package,and%20annotation%20are%20also%20supported).
2. Please install R package "tibble". More information can be found [here](https://tibble.tidyverse.org/)
3. Please install R package "org.Hs.eg.db". More information can be found [here](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
4. Please install R package "TxDb.Hsapiens.UCSC.hg38.knownGene". More information can be found [here](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html)
### Python Environment Set up
Python environment can be easily set up with [Anaconda](https://www.anaconda.com/). After `Ananconda` is installed in your OS, create a new environment from the dependencies recorded in the `environment.yml` file in this repository:
```bash
# under the root folder of the repository
conda create -f environment.yml -n mtlrank_env
```
The `mtlrank_env` is the new environment name for hosting all dependencies. Replace it with any name you prefer. Use the following command to activate the environment each time before running any scripts in this repository.
```bash
conda activate mtlrank_env
```
## Script files
- **preprocess_tfActivity_average.ipynb**: Preprocessing script for computing TF activity scores.
- **preprocess_tfExpressions.ipynb**: Preprocessing script for computing TF RPKM expressions.
- **preprocess_tfMarker_data.ipynb**: Preprocessing script for data downloaded from TF-Marker database.
- **CV_MTLRANK_expTFAPredVelo.ipynb**: **cross-validation** training and testing for MTLRank method using TF expressions + TF activities as inputs and target gene velocity values as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **CV_MTLRANK_expTFAPredVelo_SNARE.ipynb**: **cross-validation** training and testing for MTLRank method using TF expressions + TF activities as inputs and target gene velocity values as outputs. Data was generated from SNARE-seq and ChIP-seq.
- **CV_MTLRANK_expTFAPredVelo_SNAREAverage.ipynb**: **cross-validation** training and testing for MTLRank method using TF expressions + TF activities as inputs and target gene velocity values as outputs. Data was generated from SNARE-seq and ChIP-seq. TF activites computed from SNARE-seq data was averaged across cells.
- **CV_basenn_expTFAPredExp.ipynb**: **cross-validation** training and testing for baseline neural network method using TF expressions + TF activities as inputs and target gene expressions as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **CV_basenn_expTFAPredVelo.ipynb**: **cross-validation** training and testing for baseline neural network method using TF expressions + TF activities as inputs and target gene velocity values as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **CV_genie3_expTFAPredExp.ipynb**: **cross-validation** training and testing for GENIE3 method using TF expressions + TF activities as inputs and target gene expressions as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **CV_grnboost2_expTFAPredExp.ipynb**: **cross-validation** training and testing for Grnboost2 method using TF expressions + TF activities as inputs and target gene expressions as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **CV_lasso_expTFAPredVelo.ipynb**: **cross-validation** training and testing for Grnboost2 method using TF expressions + TF activities as inputs and target gene velocity values as outputs. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **MTLRANK_buildModels.ipynb**: **Model training** that uses all available data for MTLrank framework. Data was generated from scRNA-seq, scATAC-seq, and ChIP-seq.
- **MTLRANK_computeShap.ipynb**: **TF ranking based on trained models**. Ranking was performed for each individual model.
- **annotate_peaks.R**: R script for annotating ChIP-seq peaks.
- **generate_network_validate_TF.ipynb**: Script for **1)**generating the final tissue-specific networks from TF ranking results and **2)**validating the identified TFs by data from [TF-Marker](http://bio.liclab.net/TF-Marker/) database.
