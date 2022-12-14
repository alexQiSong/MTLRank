import os
import re
import subprocess

def download():
    # Pull singularity image for chipseeker
    print("Downloading singularity image for ChIPseeker..")
    if not os.path.exists("singularity"):
        os.mkdir("singularity")
    commands = "singularity pull singularity/ubuntu-18.04-R3.6.3-ChIPseeker.sif library://alexsong0374/chipseeker/chipseeker:latest"
    subprocess.run(commands, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Done")
    # File URLs.
    scrna_urls = {"Liver":[
                        "https://assets.hubmapconsortium.org/7cc96210e80d5736ca38711503e4ff06/scvelo_annotated.h5ad",
                        "https://assets.hubmapconsortium.org/7cc96210e80d5736ca38711503e4ff06/secondary_analysis.h5ad",
                        "https://assets.hubmapconsortium.org/f7afbf7de3dc1367582eebb27ae130b7/scvelo_annotated.h5ad",
                        "https://assets.hubmapconsortium.org/f7afbf7de3dc1367582eebb27ae130b7/secondary_analysis.h5ad",
                        "https://assets.hubmapconsortium.org/0e184744c4443f859984bbd3b7e95d24/scvelo_annotated.h5ad",
                        "https://assets.hubmapconsortium.org/0e184744c4443f859984bbd3b7e95d24/secondary_analysis.h5ad",
                        "https://assets.hubmapconsortium.org/bcea4743ccb03d1c4f51c82b88a528d6/scvelo_annotated.h5ad",
                        "https://assets.hubmapconsortium.org/bcea4743ccb03d1c4f51c82b88a528d6/secondary_analysis.h5ad"
                    ],
                "Spleen":[
                    "https://assets.hubmapconsortium.org/319aafd9420a0c1c6f175d6a2ef060a9/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/319aafd9420a0c1c6f175d6a2ef060a9/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/29a538c3ddb396dee26188ae1151da46/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/29a538c3ddb396dee26188ae1151da46/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/f412e76986c1012ea9589d545ed8f043/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/f412e76986c1012ea9589d545ed8f043/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/9a36e5319429ec6aca5a8a9fef401929/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/9a36e5319429ec6aca5a8a9fef401929/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/d48b7990d638dbede870ae9c1976e475/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/d48b7990d638dbede870ae9c1976e475/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/046251c94ea0e79ee935dd3de57e093c/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/046251c94ea0e79ee935dd3de57e093c/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/d9780d3f4eb9edfe275abaa32ff8633b/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/d9780d3f4eb9edfe275abaa32ff8633b/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/7b99b447ffc977a3f6f890d32c7238b3/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/7b99b447ffc977a3f6f890d32c7238b3/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/d9780d3f4eb9edfe275abaa32ff8633b/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/d9780d3f4eb9edfe275abaa32ff8633b/scvelo_annotated.h5ad",
                    "https://assets.hubmapconsortium.org/391061b538480d6d630bdfec283c0293/secondary_analysis.h5ad",
                    "https://assets.hubmapconsortium.org/391061b538480d6d630bdfec283c0293/scvelo_annotated.h5ad"
                ]
           }

    scatac_urls = {"Liver":[   
                            "https://assets.hubmapconsortium.org/73471388c8fb65d21f964a1df408db1f/peaks.combined.bed",
                            "https://assets.hubmapconsortium.org/73471388c8fb65d21f964a1df408db1f/cell_by_gene.hdf5",
                            "https://assets.hubmapconsortium.org/7e6657cc565a7d89a9041dece3816ffa/peaks.combined.bed",
                            "https://assets.hubmapconsortium.org/7e6657cc565a7d89a9041dece3816ffa/cell_by_gene.hdf5"
                        ],
                   "Spleen":[
                             "https://assets.hubmapconsortium.org/880300ed65bc837c3a81aeea474e395c/peaks.combined.bed",
                             "https://assets.hubmapconsortium.org/880300ed65bc837c3a81aeea474e395c/cell_by_gene.hdf5"
                        ]
        }

    tf_list_url = "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"

    # Create data folders
    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists(f"data/rna"):
        os.mkdir(f"data/rna")
    if not os.path.exists(f"data/atac"):
        os.mkdir(f"data/atac")
    if not os.path.exists(f"data/chipseq_bed"):
        os.mkdir(f"data/chipseq_bed")
    if not os.path.exists(f"data/chipseq_qc"):
        os.mkdir(f"data/chipseq_qc")
    if not os.path.exists(f"data/tf_list"):
        os.mkdir(f"data/tf_list")
    if not os.path.exists(f"data/gtf"):
        os.mkdir(f"data/gtf")

    # Download tissue-specific scRNA-seq and scATAC-seq data
    print("===========================")
    for tissue in scrna_urls.keys():

        if not os.path.exists(f"data/rna/{tissue}"):
            os.mkdir(f"data/rna/{tissue}")
        if not os.path.exists(f"data/atac/{tissue}"):
            os.mkdir(f"data/atac/{tissue}")

        print(f"Downloading data for tissue: {tissue}...")

        for scrna_url in scrna_urls[tissue]:
            ID = re.match(".*/(.*)/[^/]+$",scrna_url).group(1)
            filename = re.match(".*/.*/([^/]+)$",scrna_url).group(1)

            if not os.path.exists(f"data/rna/{tissue}/{ID}"):
                os.mkdir(f"data/rna/{tissue}/{ID}")

            cmd = f"wget -O data/rna/{tissue}/{ID}/{filename} {scrna_url}"
            res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"{filename}(scRNA-seq) in {ID} finished")

        for scatac_url in scatac_urls[tissue]:
            ID = re.match(".*/(.*)/[^/]+$",scatac_url).group(1)
            filename = re.match(".*/.*/([^/]+)$",scatac_url).group(1)

            if not os.path.exists(f"data/atac/{tissue}/{ID}"):
                os.mkdir(f"data/atac/{tissue}/{ID}")

            cmd = f"wget -O data/atac/{tissue}/{ID}/{filename} {scatac_url}"
            subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"{filename}(scATAC-seq) in {ID} finished")

        print(f"Tissue: {tissue} finished.")
        print("===========================")

    # Download TF list
    print(f"Downloading human TF list...")
    cmd = f"wget -O data/tf_list/tf_list.csv {tf_list_url}"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(f"Human TF list download finished")
    print("===========================")

    # Download GTF file
    print(f"Downloading human genome GTF file...")
    cmd = f"wget -O data/gtf/gencode.v32.annotation.gtf.gz {gtf_url} && gzip -d data/gtf/gencode.v32.annotation.gtf.gz"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(f"Human GTF file download finished")
    print("===========================")

    note = '''
    NOTE: Due to data sharing policy, we are not allowed to directly share ChIP-seq data from CistromeDB.
    Please go to http://cistrome.org/db/#/bdown and select Human_Factor to download all bed files
    and a QC file for ChIP-seq data. All bed files are downloaded as a gz file named "human_factor.tar.gz".
    Download and decompress this file and move all bed files to data/chipseq_bed/. QC file is named as
    "human_factor_full_QC.txt". Download this file and move it to data/chipseq_qc/
    '''
    print(note)
    return None

if __name__ == "__main__":
    download()