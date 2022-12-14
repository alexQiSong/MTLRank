# Load modules
import pandas as pd
import numpy as np
import os
import h5py
import anndata

from collections import defaultdict
from joblib import Parallel, delayed

def get_expressions():
    #------------------------------------------------------------------------------------------
    # Step1. Load genome annotation file (for getting gene length information to compute rpkm)
    #------------------------------------------------------------------------------------------
    anno = pd.read_csv("data/gtf/gencode.v32.annotation.gtf",
                       sep="\t",
                       comment = "#",
                      header = None)
    anno = anno.loc[anno[2] == "gene"]
    anno = pd.DataFrame(
        {
            'gene_ensembl':anno[8].str.extract("gene_id \"([^\.]+)\.")[0].values,
            'gene_length':(anno[4] - anno[3]).values
        }
    )
    anno.index = anno["gene_ensembl"]
    anno = anno.drop_duplicates()

    '''
    This shows sample IDs for scRNA-seq data. You may download additonal datasets from the HuBMap data portal: https://portal.hubmapconsortium.org/. Alternatively, you may perfrom batch download datasets using commandline tool as described in https://software.docs.hubmapconsortium.org/clt/index.html. This script assumes you have run download_data.py to download the sample data or you have formatted your own input data in the same way the sample data was processed.
    '''
    rootdir = "data/rna/"

    # Velocity matrix
    v_mat = defaultdict(list)

    # Expression matrix
    e_mat = defaultdict(list)

    tissues = ["Liver","Spleen"]

    # Iterate over tissues
    for tissue in ["Liver","Spleen"]:
        if os.path.isdir(f"{rootdir}/{tissue}"):
            print(f"Data(ID) for {tissue}:")
            for Dir in os.listdir(rootdir + tissue):
                velo_path = "/".join([rootdir,tissue,Dir,"scvelo_annotated.h5ad"])
                expr_path = "/".join([rootdir,tissue,Dir,"secondary_analysis.h5ad"])
                if os.path.exists(velo_path):
                    print(f"{Dir}")

    #-----------------------------------------------------------------------
    # Step2. Read velocity and expressions. Compute RPKMs from expressions
    #-----------------------------------------------------------------------
    rootdir = "data/rna"

    # Velocity matrix
    v_mat = defaultdict(list)

    # Expression matrix
    e_mat = defaultdict(list)

    # Iterate over tissues
    for tissue in os.listdir(rootdir):
        print(f"Processing expressions & velocities for {tissue}")
        for Dir in os.listdir(f"{rootdir}/{tissue}"):
            velo_path = f"{rootdir}/{tissue}/{Dir}/scvelo_annotated.h5ad"
            expr_path = f"{rootdir}/{tissue}/{Dir}/secondary_analysis.h5ad"

            # Read velocities and expression data
            v = anndata.read_h5ad(velo_path)
            expr = anndata.read_h5ad(expr_path)
            res = np.ravel(expr.layers['spliced'].todense().sum(axis = 1))

            cells = v.obs.index.intersection(expr.obs.index).values

            # Get gene ids used for velocities
            v_genes = v.var.index[v.var["velocity_genes"] == True].values 
            v_genes = [gene.split('.')[0] for gene in v_genes] # Strip of the '.XX' for ensmebl ids

            # Filter gene ids used for expressions. Keep only the genes that have length info available
            expr.var.index = pd.Index([gene.split('.')[0] for gene in expr.var.index]) # Strip of the '.XX' for ensmebl ids

            select = expr.var.index.isin(anno.index)
            rpkm_expr = expr.layers['spliced'][:,select]
            e_genes = expr.var.index[select]

            # normalize expression data to rpkm values
            lib_size = rpkm_expr.sum(axis = 1)
            gene_length = anno.loc[e_genes,]["gene_length"].values.reshape(1,-1)
            rpkm_expr = (rpkm_expr * 1000000 * 1000) / (lib_size * gene_length) 

            # For expression matrix, keep overlapping cells in both matrices
            # For velocity matrix, keep velocity genes and overlapping cells in both matrices 
            v_select = v.layers['velocity'][v.obs.index.isin(cells),:][:,v.var["velocity_genes"] == True]
            e_select = rpkm_expr[expr.obs.index.isin(cells),:]


            v_mat[tissue].append(pd.DataFrame(v_select,
                                            index = cells,
                                            columns = v_genes
                                            )
                                )

            e_mat[tissue].append(pd.DataFrame(e_select,
                                            index = cells,
                                            columns = e_genes
                                            )
                                )
        v_mat[tissue] = pd.concat(v_mat[tissue])
        e_mat[tissue] = pd.concat(e_mat[tissue], join = 'inner')
        print(f"Processing for {tissue} finished")

    #----------------------------------------
    # Step3. Save velocity and RPKM matrices
    #----------------------------------------
    if not os.path.exists("data/processed/"):
        os.mkdir("data/processed/")
    if not os.path.exists("data/processed/"):
        os.mkdir("data/processed/")

    with h5py.File("data/processed/rpkm.hdf5", "w") as f:
        for tissue,mat in e_mat.items():
            f.create_dataset("{}/exp".format(tissue), data = mat.values)
            f.create_dataset("{}/ensembl".format(tissue), data = mat.columns.values.astype('S'))
            f.create_dataset("{}/barcode".format(tissue), data = mat.index.values.astype('S'))

    with h5py.File("data/processed/velo.hdf5", "w") as f:
        for tissue,mat in v_mat.items():

            f.create_dataset("{}/velo".format(tissue), data = mat.values)
            f.create_dataset("{}/ensembl".format(tissue), data = mat.columns.values.astype('S'))
            f.create_dataset("{}/barcode".format(tissue), data = mat.index.values.astype('S'))
    return None

