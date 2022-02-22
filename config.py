import os

app_dictionary = {
    "fastqc": "quay.io/biocontainers/fastqc@sha256:319b8d4eca0fc0367d192941f221f7fcd29a6b96996c63cbf8931dbb66e53348",
    "multiqc": "quay.io/biocontainers/multiqc@sha256:82dae6463e1b19fafb6022401186300b66decf5ce319a725271700fe4e32e12a",
    "pycoqc": "quay.io/biocontainers/pycoqc@sha256:ea0a084751a0b48b5ffe90e9d3adfa8f57473709a1b0a95c9cb38d434ee3a9a2",
    "minimap2":
        "quay.io/biocontainers/minimap2@sha256:7f95eecc8eeee8ef8ae7e24d1d1a49ee056fb21d72aea4e2def97484f8a206c5",
    "samtools":
        "quay.io/biocontainers/samtools@sha256:2b911396b907769945a5446c8779d7be83f999bb4211733a8259040e67f4065f",
    "flye": "quay.io/biocontainers/flye@sha256:f895c72298ea3ae568c265cfb575fefeca768c42870cfea0ef3a4cfc35704086",
    "mlst": "quay.io/biocontainers/mlst@sha256:d8d8e731f165df398d95ad0969333aef45bcae678da11a26a1d2e76d44bc2698",
    "centrifuge":
        "quay.io/biocontainers/centrifuge@sha256:e3ce6d3d83a1df5327ee27b66d4c9eedb05b7bcd2eae8f78a7f7b9c1e8672c1c",
    "pyfaidx":
        "quay.io/biocontainers/pyfaidx@sha256:96bfff4ed96cb9149c54d21d6ea8d6d9eadcb4c0838b8e27df4b929cb5fe4e2f"
}

refseq_dict = {
    "Acinetobacter_baumannii": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Acinetobacter_baumannii/representative"
                               "/GCF_008632635.1_ASM863263v1/GCF_008632635.1_ASM863263v1_genomic.fna.gz",
    "Enterobacter_cloacae_complex":
        "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterobacter_cloacae/representative/"
        "GCF_000770155.1_ASM77015v1/GCF_000770155.1_ASM77015v1_genomic.fna.gz",
    "Escherichia_coli": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/"
                        "GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz",
    "Klebsiella_pneumoniae": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Klebsiella_pneumoniae/reference/"
                             "GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz",
    "Pseudomonas_aeruginosa": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Pseudomonas_aeruginosa/reference/"
                              "GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
    "Serratia_marcescens": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Serratia_marcescens/representative/"
                           "GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_genomic.fna.gz",
    "Staphylococcus_aureus": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/"
                             "GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
    "Enterococcus_faecium": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterococcus_faecium/representative/"
                            "GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.fna.gz"
}