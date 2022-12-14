import argparse
import pandas as pd
import numpy as np
import os

def get_tf_activities(n_jobs):
    #----------------------------------------------------------------------------
    # Step 1. Read the QC report file for all TFs and filter out low quality TFs.

    # You will need to download the QC table (named as "human_factor_full_QC") by yourself from CistromeDB: http://cistrome.org/db/#/bdown
    qc = pd.read_csv("data/chipseq_qc/human_factor_full_QC.txt",sep = "\t")

    '''
    Transcription factor selection based on quality thresholds.
    Please refer to http://cistrome.org/chilin/_downloads/instructions.pdf for more info.
    1. FastQC is the sample’s median sequence quality scores. ChiLin calculates these scores using the FastQC software[2]. A good sequence quality score is ≥ 25.
    2. Original total reads is the sample’s raw total reads number.
    3. Uniquely mapped reads is the number of reads with mapping quality above 1. First, ChiLin aligns reads onto user-specified genomes. Then, it filters the SAM files. The uniquely mapped RATIO is the uniquely mapped reads divided by the total reads. **A good uniquely mapped ratio is ≥ 60%.**
    4. Unique locations of 4M reads is the number of genomic locations with one or more uniquely mapped reads (unique locations) from sub-sampled 4M reads. Unique locations ratio unique locations number divided by total number of uniquely mapped reads. ChiLin estimates NRF by dividing the number of unique locations by 4M sampled uniquely mapped reads. If reads are less than 4M, then ChiLin uses the total reads instead. ChiLin reports number of unique locations and the unique locations ratio. A good unique locations of 4M reads should be ≥ 70%.
    5. Locations with only 1 read from 4M reads number (ratio) is the number of locations with read number equal to 1 (N1). The ratio is N1 divided by 4M reads unless the total reads is less than 4M, in which case the total reads is used. A good score for this metric is > 70%.
    6. PBC of 4M reads is N1 (see 5) divided by unique locations (see 4). **A good PBC score is ≥ 80%.**
    7. Fragment size of 4M reads is in silico estimation of your size selection through maximum cross correlation. The estimation should to be close to the size selected in your experiment.
    8. Exon/DHS/Promoter ratio of 4M reads is the estimated ratio of reads falling in these regions (from a 4M reads sub-sample). Exons regions are defined as the merged exons regions from the RefSeq gene table. Promoter regions are defined as the RefSeq TSS +/- 2kb regions. Union DHS regions are called from ENCODE II UW DNase-seq Hypersensitive regions. The IP group samples should have higher reads ratios than the control group samples.
    9. **FRiP of 4M non-chrM reads is used for evaluating the signal to noise ratio.** First, ChiLin removes chrM reads from the total reads. Then ChiLin sub-samples 4M of these reads. Finally, it calculates the ratio of the sub-sample which fall under the called peaks. **A good FRiP score is ≥ 1%.**
    10. Replicates total peaks are the total peaks number called by MACS2 with fixed extension size and q value cutoff. A good peaks number depends on your experiment.
    11. Replicates 10 fold confident peaks are the number of peaks called by MACS2 where the fold change is ≥ 10. **A good number is above 500**.
    '''
    # Select transcription factors with high quality, total number of selected TFs = 623
    qc = qc.loc[(qc.FastQC >= 25) & 
     (qc.UniquelyMappedRatio >= 0.6) &
     (qc.PBC > 0.8) &
     (qc.FRiP > 0.01) &
     (qc.PeaksFoldChangeAbove10 > 500), ]

    # Change SMAD2/3 -> SMAD23 to use the TF name as file name
    qc.loc[qc["Factor"] == "SMAD2/3","Factor"] = "SMAD23"

    #----------------------------------------------------------------------------
    #Step2. Check number of cells, and number of peaks for each tissue

    import h5py
    import os

    rootdir = "data/atac"

    # Iterate over tissues, read in cell by bin matrix, and peaks files. Then count number of cells and peaks.
    for tissue in os.listdir(rootdir):

        n_peaks = 0
        n_cells = 0

        # Read peaks and cell by bin matrix
        for Dir in os.listdir(f"{rootdir}/{tissue}"):
            peaks = pd.read_csv(f"{rootdir}/{tissue}/{Dir}/peaks.combined.bed",
                            sep = "\t", header = None)
            n_peaks += peaks.shape[0]
            n_cells += h5py.File(f"{rootdir}/{tissue}/{Dir}/cell_by_gene.hdf5","r")['row_names'].shape[0]

        print(f"Tissue {tissue}, {n_cells} cells and {n_peaks} peaks in total")

    #---------------------------------------------------------------------------------------------------------
    # Step 3. Read the scATAC-seq peaks, and build scATAC-seq search tree for each tissue and each chromosome.
    from intervaltree import IntervalTree,Interval
    from collections import defaultdict
    import re

    rootdir = "data/atac/"
    search_tree = defaultdict(dict)

    # Iterate over tissues and chromosomes. For each tissue, buid interval trees.
    print("=============================================")
    print("Building search trees for scATAC-seq peaks...")
    for tissue in os.listdir(rootdir):

        peaks = []
        for Dir in os.listdir(f"{rootdir}/{tissue}"):
            peaks.append(pd.read_csv(f"{rootdir}/{tissue}/{Dir}/peaks.combined.bed", sep = "\t"))
            peaks[-1].columns = ["chrom","chromStart","chromEnd"]

        # Combine peaks from the same tissue
        peaks = pd.concat(peaks)
        peaks['chrom'] = 'chr' + peaks['chrom']

        # For each chromosome, build an interval tree for search.
        for chrom,df in peaks.groupby('chrom'):

            # Remove parenthesis and spaces in the tissue name
            tissue = re.sub('\s+',"_",tissue)
            tissue = re.sub('\(',"",tissue)
            tissue = re.sub('\)',"",tissue)

            search_tree[tissue][chrom] = IntervalTree(Interval(row.chromStart,row.chromEnd) for i,row in df.iterrows())
            search_tree[tissue][chrom].merge_overlaps()
        print(f"Search tree for {tissue} is completed")


    # Step 3. For each ChIP-seq peaks of TF, compute the overlapping regions with scATAC-seq, and map the peak to nearest genes.
    # Need to install ChIPseeker R package before performing this step. You may also find our singularity image for running ChIPseeker.

    from joblib import Parallel, delayed
    from tqdm import tqdm
    import subprocess

    if not os.path.exists("results/tf_activities_bed"):
        os.mkdir("results/tf_activities_bed")

    # Create tissue subfolder for saving TFs.
    for tissue in search_tree.keys():

        if not os.path.exists(f"results/tf_activities_bed/{tissue}"):
            os.mkdir(f"results/tf_activities_bed/{tissue}")

    # Function for searching overlapping regions between scATAC-seq peaks and ChIP-seq peaks 
    def search_overlap(tf,df):

        # Read all peak files for a TF
        peaks = []
        for idx,row in df.iterrows():
            peaks.append(
                            pd.read_csv("data/chipseq_bed/" + str(row.DCid) + "_sort_peaks.narrowPeak.bed",
                            sep = "\t",
                            header = None)
                        )
            peaks[-1].columns = ["chrom",
                             "chromStart",
                             "chromEnd",
                             "name",
                             "score",
                             "strand",
                             "signalValue",
                             "pValue",
                             "qValue",
                             "peak"]
            peaks[-1]["tf"] =  np.repeat(row.Factor,peaks[-1].shape[0])
            peaks[-1]["DCid"] = np.repeat(row.DCid,peaks[-1].shape[0])

        peaks = pd.concat(peaks)

        # Search the overlapped between scATAC-seq peaks and ChIP-seq peaks for each TF in each tissue
        # each row of 'peaks' represents a peak in ChIP-seq data.
        # each 'peaks' table represents all peaks of a TF from ChIP-seq data.
        # 1. iterate over all tissues
        for tissue,tree in search_tree.items():

            # 2. iterate over all peaks of a TF
            res = []
            for idx,row in peaks.iterrows():
                chrom = row['chrom']
                start = row['chromStart']
                end = row['chromEnd']
                sigVal = row['signalValue']

                if chrom in tree:
                    intervals = tree[chrom][start:end]
                    intervals = sorted(list(intervals), key = lambda x: x[0])
                    over_len = 0
                    for interval in intervals:
                        over_len += min(end,interval[1]) - max(start,interval[0])
                    res.append(over_len/(int(end) - int(start)))
                else:
                    res.append(None)

            # After all peaks of a TF is done. All peaks of the TF is saved as one file, specific to each tissue.
            peaks["chip_atac_weight"] = res
            outname = f"results/tf_activities_bed/{tissue}/{tf}.csv".format(tissue,tf)
            peaks.to_csv(outname)

            # map TF peaks to nearest genes by R ChIPseeker
            commands = f"singularity exec singularity/ubuntu-18.04-R3.6.3-ChIPseeker.sif Rscript annotate_peaks.R {outname} {outname}"
            subprocess.run(commands, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        return(None)

    print("======================================================================")
    print("Searching overlapping regions for each TF. This may take some time...")
    res = Parallel(n_jobs=n_jobs)(delayed(search_overlap)(tf,df) 
                                  for i,(tf,df) in enumerate(qc.groupby("Factor")))
    print("Search is completed")


    # Step 4. For each tissue, read the results and generate the TF-gene activity matrix
    # Get non-zero and non-na TF-gene pairs

    import os
    from tqdm import tqdm
    from collections import defaultdict

    folder = "results/tf_activities_bed"
    tissues = os.listdir(folder)
    tf_gene = defaultdict(list)

    print("======================================================================")
    print("Generating tf-gene pairs...")

    for tissue in tissues:
        tf_files = os.listdir(f"{folder}/{tissue}")

        print(f"Generating tf_gene pairs for {tissue}...")
        for i in tqdm(range(len(tf_files))):
            df = pd.read_csv(f"{folder}/{tissue}/{tf_files[i]}")
            tf_gene[tissue].append(df.loc[~df["gene"].isna() & (df["chip_atac_weight"] != 0),["tf","gene","chip_atac_weight"]])
        tf_gene[tissue] = pd.concat(tf_gene[tissue])
    print("Done")


    # Get mean of weights for each TF-gene pair
    print("=====================================================")
    tf_gene_mean = defaultdict(list)
    for tissue,tab in tf_gene.items():
        print(f"Computing TF activity mean scores for {tissue}...")
        tf_gene_mean[tissue] = tab.groupby(["tf","gene"], as_index = False)['chip_atac_weight'].mean()
    print("Done")


    # Get sum of weights for each TF-gene pair
    print("=====================================================")
    tf_gene_sum = defaultdict(list)
    for tissue,tab in tf_gene.items():
        print(f"Computing TF activity sum scores for {tissue}...")
        tf_gene_sum[tissue] = tab.groupby(["tf","gene"], as_index = False)['chip_atac_weight'].sum()
    print("Done")


    # Convert data frames to tf-gene activities matrices
    print("=====================================================")
    print("Saving TF activity scores...")
    tf_gene_mat_mean = dict()
    for tissue,mat in tf_gene_mean.items():
        tf_gene_mat_mean[tissue] = mat.pivot(index = "tf", columns = "gene", values = "chip_atac_weight")

    tf_gene_mat_sum = dict()
    for tissue,mat in tf_gene_sum.items():
        tf_gene_mat_sum[tissue] = mat.pivot(index = "tf", columns = "gene", values = "chip_atac_weight")


    # Save tf-gene activities matrix
    if not os.path.exists("results/tf_activities_mat/"):
        os.mkdir("results/tf_activities_mat/")

    for tissue in tf_gene_mat_mean.keys():
        tf_gene_mat_mean[tissue].to_csv(f"results/tf_activities_mat/{tissue}_tfGeneMean.csv")

    for tissue in tf_gene_mat_sum.keys():
        tf_gene_mat_sum[tissue].to_csv(f"results/tf_activities_mat/{tissue}_tfGeneSum.csv")
    print("Done")
    return None