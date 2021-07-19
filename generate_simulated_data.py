import install_containers
import wget
import os
import gzip
import shutil
import glob

app_dictionary = {
    "nanosim": "quay.io/biocontainers/nanosim@sha256:d99389f4fafd8a36547cf5c2a6996a97d929482090682b1a4d070c28069d199b",
    "nanosim_h": "quay.io/biocontainers/nanosim-h@sha256:"
                 "76e3d6ab85a917623886d04b49504f1c0865dcfb6fa27cf9d8bd1a7145a26150"
}

reference_list = {"A_baumannii": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Acinetobacter_baumannii/"
                                "representative/GCF_002116925.1_ASM211692v1/GCF_002116925.1_ASM211692v1_genomic.fna.gz",
                 "E_cloacae_complex": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterobacter_cloacae_complex_sp"
                                      "./latest_assembly_versions/GCF_900322725.1_C45/GCF_900322725.1_C45_genomic.fna."
                                      "gz",
                 "E_coli": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845."
                           "2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
                 "K_pneumoniae": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Klebsiella_pneumoniae/reference/"
                                 "GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz",
                 "P_aeruginosa": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Pseudomonas_aeruginosa/reference/"
                                 "GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz",
                 "S_marcescens": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Serratia_marcescens/representative/"
                                 "GCF_001034395.1_Serr_marc_UCI87_V1/GCF_001034395.1_Serr_marc_UCI87_V1_genomic.fna.gz",
                 "S_aureus": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Staphylococcus_aureus/reference/"
                             "GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz",
                 "E_faecium": "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Enterococcus_faecium/representative/"
                              "GCF_000336405.1_ASM33640v1/GCF_000336405.1_ASM33640v1_genomic.fna.gz"
                 }

def download_refs(out_dir):
    """
    If reference sequence directory is empty, download microbial reference genomes from NCBI using reference_list
    dictionary.
    """
    if not os.listdir(out_dir):
        for key in reference_list:
            file_name = str(reference_list[key]).rsplit("/", 1)[1]
            if os.path.isfile(out_dir + file_name):
                print(key + " reference already downloaded from NCBI")
            else:
                print("\nDownloading reference sequence from NCBI for " + key)
                wget.download(reference_list[key], out=out_dir)

#def nanosim():

def unzip(refs_dir):
    """
    Unzip microbial reference genomes in specified directory.
    """
    for file in glob.glob(refs_dir + "*.gz"):
        unzipped_name = str(file).rsplit(".gz", 1)[0]
        if os.path.isfile(unzipped_name):
            print(unzipped_name + " already exists.")
        else:
            print("Unzipping " + file)
            with gzip.open(file, 'rb') as infile:
                with open(unzipped_name, 'wb') as outfile:
                    shutil.copyfileobj(infile, outfile)
            os.remove(file)

def nanosim_h():
    for reference in reference_list:
        "nanosim-h {} -p ecoli_R9_1D".format(reference)

def main():
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])
    microbial_refs_dir = "/media/data2/share/outbreak_pipeline/rduffin/outbreak_pipeline/test_data/input/" \
                         "microbial_refs/"
    download_refs(out_dir = microbial_refs_dir)
    unzip(refs_dir = microbial_refs_dir)


if __name__ == '__main__':
    main()
