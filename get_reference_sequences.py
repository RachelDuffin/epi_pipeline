"""
Script to download required reference sequences from NCBI refseq for all scripts.
"""
import requests
import os

refseq_dict = {
    "Acinetobacter_baumannii": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Acinetobacter_baumannii/representative/GCF_008632635.1_ASM863263v1/GCF_008632635.1_ASM863263v1_genomic.fna.gz",
    "Enterobacter_cloacae_complex": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterobacter_cloacae/representative/GCF_000770155.1_ASM77015v1/GCF_000770155.1_ASM77015v1_genomic.fna.gz",
    "Escherichia_coli": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz",
    "Klebsiella_pneumoniae": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Klebsiella_pneumoniae/reference/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz",
    "Pseudomonas_aeruginosa": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Pseudomonas_aeruginosa/reference/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
    "Serratia_marcescens": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Serratia_marcescens/representative/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_genomic.fna.gz",
    "Staphylococcus_aureus": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
    "Enterococcus_faecium": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterococcus_faecium/representative/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.fna.gz"
}

refseq_test_dict = {
    "Acinetobacter_baumannii": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Acinetobacter_baumannii/representative/GCF_008632635.1_ASM863263v1/GCF_008632635.1_ASM863263v1_genomic.fna.gz",
    "Enterobacter_cloacae_complex": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterobacter_cloacae/representative/GCF_000770155.1_ASM77015v1/GCF_000770155.1_ASM77015v1_genomic.fna.gz",
    "Escherichia_coli": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz",
    "Klebsiella_pneumoniae": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Klebsiella_pneumoniae/reference/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz",
    "Pseudomonas_aeruginosa": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Pseudomonas_aeruginosa/reference/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
    "Serratia_marcescens": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Serratia_marcescens/representative/GCF_003516165.1_ASM351616v1/GCF_003516165.1_ASM351616v1_genomic.fna.gz",
    "Staphylococcus_aureus": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
    "Enterococcus_faecium": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterococcus_faecium/representative/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.fna.gz",
    "Listeria_monocytogenes": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Listeria_monocytogenes/reference/GCF_000196035.1_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna.gz",
    "Bacillus_subtilis": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Bacillus_subtilis/reference/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz",
    "Saccharomyces_cerevisiae": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz",
    "Salmonella_enterica": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Salmonella_enterica/reference/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz",
    "Lactobacillus_fermentum": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/538/825/GCF_006538825.1_ASM653882v1/GCF_006538825.1_ASM653882v1_genomic.fna.gz",
    "Enterococcus_faecalis": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterococcus_faecalis/representative/GCF_000393015.1_Ente_faec_T5_V1/GCF_000393015.1_Ente_faec_T5_V1_genomic.fna.gz",
    "Cryptococcus_neoformans": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/045/GCF_000091045.1_ASM9104v1/GCF_000091045.1_ASM9104v1_genomic.fna.gz",
}

def download_sequences(dictionary, out_dir):
    """
    Download all sequences in supplied dictionary.

    Checks if download already exists. If it doesn't, download the file, get the file name by splitting the dictionary
    value, and append bacterial name to start using the dictionary key.
    """
    base_path = os.getcwd()

    for key in dictionary:
        file_name = str(dictionary[key]).rsplit("/", 1)[1]
        filepath= "{}/{}{}_{}".format(base_path, out_dir, key, file_name)

        if not os.path.exists(filepath):
            print("Downloading: {} refseq file".format(key))
            request = requests.get(dictionary[key])
            open(filepath, 'wb').write(request.content)
            print("Download complete: {} refseq file".format(key))

def main():
    out_dir = "/input/reference_sequences/"
    download_sequences(refseq_test_dict, out_dir)
    download_sequences(refseq_dict, out_dir)

if __name__ == '__main__':
    main()
