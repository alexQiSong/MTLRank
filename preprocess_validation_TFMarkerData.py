import pandas as pd
import numpy as np

# read ensembl to symbol mapping
ensembl_to_symbol = pd.read_csv("data/raw/id_mapping/ensembl_to_symbol.csv",index_col = 0)
ensembl_to_symbol = ensembl_to_symbol.loc[~ensembl_to_symbol["ensembl_id"].duplicated(),:]
ensembl_to_symbol.index = ensembl_to_symbol["ensembl_id"]

# Read TFMarker file, downloaded from http://bio.liclab.net/TF-Marker/
table = pd.read_csv("data/raw/tfmarker/tfMarker.txt", sep = "\t")

# Convert TF markers and their interacting genes to ensembl IDs
new_tab = []
for idx,row in table.iterrows():
    source = row["Gene Name"]
    if row["Interacting Gene"] is not np.nan:
        targets = row["Interacting Gene"].split(";")
    else:
        targets = ["NA"]
    tissue = row["Tissue Type"]
    exp_met = row["Experimental Method"]
    pmid = row["PMID"]
    for tar in targets:
        new_tab.append([source, tar, tissue, exp_met,pmid])
new_tab = pd.DataFrame(new_tab, columns = ["source","target","tissue","experiment","PMID"])
source_ensembl = new_tab.merge(ensembl_to_symbol,
                              how = 'left',
                              left_on = 'source',
                             right_on = 'gene_symbol')["ensembl_id"]
target_ensembl = new_tab.merge(ensembl_to_symbol,
                              how = 'left',
                              left_on = 'target',
                             right_on = 'gene_symbol')["ensembl_id"]

# Keep only the genes that can be mapped to ensmebl ids
select = (~pd.isna(source_ensembl)) & (~pd.isna(target_ensembl))

new_tab["source_ensembl"] = source_ensembl
new_tab["target_ensembl"] = target_ensembl
new_tab = new_tab.loc[select,]

# Save converted TF markers and their converted targets.
new_tab.to_csv("data/raw/tfmarker/tfMarker_ensmebl.csv", index = None)