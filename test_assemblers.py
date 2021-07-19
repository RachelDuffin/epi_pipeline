import install_containers
import subprocess
import os

app_dictionary = {
    "flye": "quay.io/biocontainers/flye@sha256:f895c72298ea3ae568c265cfb575fefeca768c42870cfea0ef3a4cfc35704086",
    "canu": "quay.io/biocontainers/canu@sha256:b48b52afc355477015ef60bebded3b4ab3d3099bbf5698879de8eb600c9ff1a4",
    "raven":
        "quay.io/biocontainers/raven-assembler@sha256:3bc4cc61483cc48243f6b416eaae41f24eb95f75b7a2770f8062c75b5ac53da3",
    "minimap2":
        "quay.io/biocontainers/minimap2@sha256:7f95eecc8eeee8ef8ae7e24d1d1a49ee056fb21d72aea4e2def97484f8a206c5",
    "racon": "quay.io/biocontainers/racon@sha256:16b6aeb33ee2ac726d313a9de3788adff305e24c07bfeee9f1800c060627b414",
    "medaka": "quay.io/biocontainers/medaka@sha256:6aa52d718af0f48cf6630e88b22cd7187770bbf028ef89bc54ec7fad2ff7a35f"
}

base_path = "/media/data2/share/outbreak_pipeline/rduffin/outbreak_pipeline/test_data/"
input_filepaths = {
    "input/mock_microbial_community/Zymo-GridION-EVEN-BB-SN.fq.gz": "output/assemblers/mock_microbial_community/",
    "input/mock_microbial_community/Zymo-GridION-LOG-BB-SN.fq.gz": "output/assemblers/mock_microbial_community/",
    "input/enterococcus_faecium/enterococcus/ef1_bc_75/210612_EF_R1_barcode01.fq":
        "output/assemblers/monomicrobial_samples/",
    "input/enterococcus_faecium/enterococcus/ef1_bc_75/210612_EF_R1_barcode02.fq":
        "output/assemblers/monomicrobial_samples/",
    "input/enterococcus_faecium/enterococcus/ef1_bc_75/210612_EF_R1_barcode03.fq":
        "output/assemblers/monomicrobial_samples/",
    "input/enterococcus_faecium/enterococcus/ef1_bc_75/210612_EF_R1_barcode04.fq":
        "output/assemblers/monomicrobial_samples/",
    "input/enterococcus_faecium/enterococcus/ef1_bc_75/210612_EF_R1_barcode05.fq":
        "output/assemblers/monomicrobial_samples/"
}

# Script to compare performance, scalability and accuracy between assembly tools

# PERFORMANCE AND SCALABILITY:
# Comparison of performance and memory usage between tools (single thread), using /usr/bin/time -v. Compare by:
# - Wall clock time - elapsed time from start to end of program
# - CPU time - Time the CPU spends in user mode (running program's code) and kernel mode (executing system calls)
# - Peak memory usage - max memory used by program during lifetime
# - Parallel speedup - ratio of time to run program with one thread to time to run with N threads
# This is run multiple times changing the number of threads and observing corresponding changes in metrics
# (all tools - flye, canu, raven - have multithreading support). The linux box has 20 cores.

# ACCURACY:
# Comparing the assembly to the reference genome, using the MUMner package. Compare by:
# - No. bases in assembly (should be length of ref. genome)
# - No. contigs (should be 1)
# - Average identity - % similarity between assembly and reference genome (should be near to 100%)
# - Coverage - ratio of no. aligned bases in ref to length of ref (should be near to 100%)
# - No. mismatches - no. single-base differences between assembly and ref (should be 0)
# - No. indels - number of insertions and deletions between assembly and ref

def mummer():
    """
    Accuracy comparison using MUMmer package using monomicrobial Nano-Sim H generated reads.
    """

def flye_assembly(input_filepath, sample_name, out_dir, threads):
    """
    Run test datasets through flye with differing numbers of threads to assess scalability, performance and accuracy.
    """
    # scalability/performance
    print("--------------------------\nFLYE DE NOVO ASSEMBLY for {}\n--------------------------".format(input_filepath))

    # overwrite old file
    out_dir = out_dir + "flye/"
    time_file = out_dir + "{}_flye_{}_thread.txt".format(sample_name, threads)
    time_output = open(time_file, 'w')
    time_output.close()
    with open(time_file, 'a') as filetowrite:
        stdout, stderr = subprocess.Popen("/usr/bin/time -v sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} flye "
                                          "--nano-raw {} --out-dir {} --meta --threads "
                                          "{}".format(app_dictionary["flye"], input_filepath, out_dir, threads),
                                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
        filetowrite.write(stdout.decode('utf-8'))
        os.rename("{}assembly.fasta".format(out_dir), "{}{}_flye_{}_thread.fasta".format(out_dir, sample_name, threads))
    print("--------------------------")


def canu_assembly(input_filepath, sample_name, out_dir, threads):
    """
    Run test datasets through canu with differing numbers of threads to assess scalability, performance and accuracy.
    """
    # scalability/performance
    print("--------------------------\nCANU DE NOVO ASSEMBLY for {}\n--------------------------".format(input_filepath))

    # overwrite old file
    out_dir = out_dir + "canu/"
    time_file = out_dir + "{}_canu_{}_thread.txt".format(sample_name, threads)
    time_output = open(time_file, 'w')
    time_output.close()
    assembly_prefix = "{}_canu_{}_thread".format(sample_name, threads)
    if 'Zymo' in input_filepath:
        genomeSize = "6.1944m"
    if 'enterococcus_faecium' in input_filepath:
        genomeSize = "2.85m"
    with open(time_file, 'a') as filetowrite:
        stdout, stderr = subprocess.Popen("/usr/bin/time -v sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} canu"
                                          " -p {} -d {} genomeSize={} minThreads={} maxThreads={} obtovlThreads={} "
                                          "utgovlThreads={} corThreads={} redThreads={} batThreads={} "
                                          "maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high "
                                          "corMinCoverage=0 redMemory=24 oeaMemory=24 "
                                          "-nanopore {}".format(app_dictionary["canu"], assembly_prefix, out_dir,
                                                                genomeSize, threads, threads, threads, threads, threads,
                                                                threads, threads, input_filepath),
                                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
        filetowrite.write(stdout.decode('utf-8'))
       # os.rename("{}assembly.fasta".format(out_dir), "{}{}_canu_{}_thread.fasta".format(out_dir, sample_name, threads))
    print("--------------------------")


def raven_assembly(input_filepath, sample_name, out_dir, threads):
    """
    Run test datasets through raven with differing numbers of threads to assess scalability, performance and accuracy.
    """
    # scalability/performance
    print("--------------------------\nRAVEN DE NOVO ASSEMBLY for {}\n--------------------------".format(input_filepath))

    # overwrite old file
    out_dir = out_dir + "raven/"
    time_file = out_dir + "{}_raven_{}_thread.txt".format(sample_name, threads)
    time_output = open(time_file, 'w')
    time_output.close()
    with open(time_file, 'a') as filetowrite:
        stdout, stderr = subprocess.Popen("/usr/bin/time -v sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} raven"
                                          " {} -t {}".format(app_dictionary["raven"], input_filepath, threads),
                                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
        filetowrite.write(stdout.decode('utf-8'))
        os.rename("{}assembly.fasta".format(out_dir), "{}{}_raven_{}_thread.fasta".format(out_dir, sample_name, threads))
    print("--------------------------")


def main():
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])

    # conduct comparison of metagenomic assemblers, using 5 E. faecium isolates (monomicrobial samples), and Zymo
    # mock community data
    for fq in input_filepaths:
        sample_name = str(fq).rsplit("/", 1)[1].rsplit(".")[0]
        input_filepath = base_path + fq
        out_dir = base_path + input_filepaths[fq]
        for threads in [8, 4, 2, 1]:
            raven_assembly(input_filepath, sample_name, out_dir, threads)
            flye_assembly(input_filepath, sample_name, out_dir, threads)
            canu_assembly(input_filepath, sample_name, out_dir, threads)

if __name__ == '__main__':
    main()
