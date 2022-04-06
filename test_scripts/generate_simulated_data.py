"""
Generates simulated monomicrobial and metagenomic datasets from a set of reference genome sequences in a reference
sequence directory.

N.B. This script was not used and is incomplete
"""

import os
import gzip
import shutil
import glob
import test_assemblers
import get_reference_sequences
import subprocess
import config

tool_dict = {
    "nanosim": "quay.io/biocontainers/nanosim@sha256:b3cdc9ffb111069b787c5b2f744f2db9ffe552acdf10ab34428bd4ffd9f179d0",
    "nanosim-h":
        "quay.io/biocontainers/nanosim-h@sha256:76e3d6ab85a917623886d04b49504f1c0865dcfb6fa27cf9d8bd1a7145a26150"
}

# median total length of each organism's genome - from refseq. 9/10/2021
median_genome_length_dict = {
    "Acinetobacter_baumannii": 3.97473,
    "Enterobacter_cloacae_complex": 5.03448,
    "Escherichia_coli": 5.12126,
    "Klebsiella_pneumoniae": 5.59628,
    "Pseudomonas_aeruginosa": 6.61337,
    "Serratia_marcescens": 5.21305,
    "Staphylococcus_aureus": 2.83686,
    "Enterococcus_faecium": 2.92265
}

def get_error_profile(infile, out_dir, refs_dir, base_path):
    """
    Generate error profile using test data.
    """
    # if error profile is not already downloaded, download it
    if not os.path.isdir("{}ecoli_R9_1D".format(out_dir)):
        error_command = "git clone https://github.com/karel-brinda/NanoSim-H.git && cd NanoSim-H && " \
                        "git reset --hard 51fc8cd && cd .. && " \
                        "cp -r NanoSim-H/nanosimh/profiles/ecoli_R9_1D/ ecoli_R9_1D && rm -r -f NanoSim-H"
        process = subprocess.Popen(error_command, shell=True, universal_newlines=True, bufsize=1,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("--------------------------")


def nanosim(refs_dir, reads_per_megabase, genome_length_dict, base_path, out_dir, error_profile_outdir):
    """
    Generate simulated genomic reads per reference genome.
    """
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # move to out dir as no out_dir option in command
    os.chdir(out_dir)

    for reference in glob.glob("{}*.fna".format(refs_dir)):
        print("--------------------------\nGENERATE SIMULATED GENOMIC READS {}\n"
              "--------------------------".format(reference))

        for microbe in median_genome_length_dict:
            if microbe in reference:
                # calculate number of reads for the genome (rounded to nearest integer as required for nanosim-h)
                read_count = round(genome_length_dict[microbe] * reads_per_megabase)
                # circular option used as all inputs are circular genomes

                #command = "module load apps/singularity && singularity exec --bind `pwd`:`pwd` " \
                #          "{}/nanosim-h.sif nanosim-h --circular {} -p 'ecoli_R9_1D' -o {} -n {} --max-len 53793 " \
                #          "--min-len 65".format(base_path, reference, microbe, read_count)

                command = "module load apps/singularity && singularity exec --bind `pwd`:`pwd` " \
                          "{}/nanosim.sif simulator.py genome -rg {} -c {} -n {} -max 53793 " \
                          "-min 65 -b guppy -dna_type circular -o {} -t 16 --fastq".format(base_path, reference,
                                                                                           error_profile_outdir,
                                                                                           read_count, microbe)
                print(command)
                run_process(command)

    os.chdir(base_path)
    print("--------------------------")

def run_process(command):
    """
    Runs command as subprocess
    """
    process = subprocess.Popen(command, shell=True, universal_newlines=True, bufsize=1,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = process.communicate()

def nanosim_metagenome(refs_dir):
    """
    Generate simulated metagenomic reads per reference genome.
    """
    out_dir = "/media/data2/share/outbreak_pipeline/rduffin/outbreak_pipeline/test_data/input/" \
              "simulated_metagenomic_reads"
    os.chdir(out_dir)
    command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} nanosim simulator.py metagenome -gl {} -n {} " \
                  "-max 53793 -min 65 -med 1503 -b guppy".format(tool_dict[nanosim], genome_list, read_number)
    # create genome list file, use the calculated max/min read lengths from my other script here?? although were they
    # calculated from enterococcus?
    run_process(commmand)
    os.chdir(cwd)

def main():

    base_path = os.getcwd()
    # reads per megabase calculated using the average read count from test data, and the median read length for
    # e.faecium from refseq (9/10/2021). 223239/2.92265 = 76382 reads per megabase
    reads_per_megabase = 76382
    for tool in tool_dict:
        test_assemblers.singularity_pull(tool, tool_dict[tool])

    # Download reference sequences
    microbial_refs_dir = "{}/input/reference_sequences/".format(base_path)
    input_dir = "/users/k1927570/outbreak_pipeline/msc_project/input/"
    out_dir = "{}/input/simulated_monomicrobial_reads/".format(base_path)
    error_profile_outdir = "{}enterococcus_error_profile".format(out_dir)
    e_faecium_dir = "{}enterococcus_faecium/enterococcus/ef1_bc_75/".format(input_dir)

    error_profile_input = "/users/k1927570/outbreak_pipeline/msc_project/input/enterococcus_faecium/enterococcus/" \
                          "ef1_bc_75/210612_EF_R1_barcode07.fq"

    get_reference_sequences.download_sequences(get_reference_sequences.refseq_dict, microbial_refs_dir, base_path)

    get_error_profile(error_profile_input, error_profile_outdir, microbial_refs_dir, base_path)

    nanosim(microbial_refs_dir, reads_per_megabase, median_genome_length_dict, base_path, out_dir, error_profile_outdir)
    #nanosim_metagenome(refs_dir = microbial_refs_dir)

if __name__ == '__main__':
    main()
