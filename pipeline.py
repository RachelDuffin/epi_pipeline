import glob
import os
import re
import subprocess
import sys

import install_containers


# SETTINGS
def get_identifier(file, string):
    """
    Gets the run id from the run's final_summary file
    """
    with open(file, "rt") as file_contents:
        lines = file_contents.readlines()
        for line in lines:
            if line.__contains__(string):
                selected = (line.split("=", 1))[1]
    if ".fq" in file:
        selected = selected.split()[-1].split("=")[-1]
    identifier = selected
    return identifier


def create_directory(parent_directory, directory):
    """
    Create output directories
    """
    if os.path.exists(parent_directory + "/" + directory):
        print(directory + " directory already exists")
    else:
        os.mkdir(os.path.join(parent_directory, directory))
        print("output folder " + directory + " created")


def fastqc(data_dir, out_dir, run_id):
    """
    FastQC analysis per run folder
    """
    print("--------------------------\nFASTQC\n--------------------------")
    output_filename = run_id + "_24hrs_fastqc.html"
    # if output file already present, do not re-analyse
    if os.path.isfile(out_dir + "/fastqc/" + output_filename):
        print("FastQC output for run " + run_id + " already exists")
    else:
        print("Creating fastqc file for run " + run_id)
        for file in glob.glob(data_dir + "/*_24hrs.fq"):
            barcode = get_identifier(file=file, string="barcode")
            subprocess.run("module load apps/singularity; singularity exec apps/fastqc.sif fastqc " + file + " -o " +
                           out_dir + "/fastqc", shell=True)
    print("--------------------------")


def pycoqc(data_dir, out_dir):
    """
    Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes() function.
    Then creates a pycoQC json report per barcode.
    """
    print("--------------------------\nPYCOQC\n--------------------------")
    # If pycoQC not yet run, split the summary sequencing files according to barcodes and run pycoQC
    if not os.listdir(out_dir + "/pycoqc"):
        print("pycoQC not yet run")
        split_barcodes(data_dir=data_dir, out_dir=out_dir)
        for file in glob.glob(out_dir + "/pycoqc/sequencing_summary_*"):
            # get barcode name from barcode_split output file names
            barcode = file.split(".", 1)[0].split("summary_", 1)[1]
            print("Creating pycoQC json report for " + file)
            subprocess.run(
                "module load apps/singularity; singularity exec apps/pycoqc.sif pycoQC -f " + file +
                " --json_outfile " + out_dir + "/pycoqc/" + barcode + "_pycoQC_output.json", shell=True)
    else:
        print("Directory not empty - pycoQC already run")
    print("--------------------------")


def split_barcodes(data_dir, out_dir):
    for file in glob.glob(data_dir + "/sequencing_summary_*"):
        print("Splitting summary sequencing file " + file + " according to barcodes")
        subprocess.run(
            "module load apps/singularity; singularity exec apps/pycoqc.sif Barcode_split --output_unclassified " +
            "--min_barcode_percent 0.0 --summary_file " + file + " --output_dir " + out_dir + "/pycoqc",
            shell=True)


def human_read_removal(data_dir, out_dir, run_id, ref):
    """
    Removes human reads from the samples by alignment to the human reference genome.
    """
    print("--------------------------\nHUMAN READ REMOVAL\n--------------------------")
    for file in glob.glob(data_dir + "/*.fq"):
        barcode = get_identifier(file=file, string="barcode")
        out_path = out_dir + "/human_read_removal/" + run_id + "_" + barcode
        if os.path.isfile(out_path + "_aligned.sam"):
            print(barcode + " already aligned")
        else:
            # align reads to human reference genome using ont-specific parameters
            print("Aligning reads to human reference genome for " + file)
            subprocess.run(
                "module load apps/singularity; singularity exec apps/minimap2.sif minimap2 -ax map-ont " + ref + " " +
                file + " > " + out_path + "_aligned.sam", shell=True)
        if os.path.isfile(out_path + "_unmapped.fastq"):
            print("Human read removal already complete for " + barcode)
        else:
            # import SAM to BAM as @SQ lines present in header, only import unassigned reads (non-human)
            print("Import non-human reads as BAM file for " + out_path + "_aligned.sam")
            subprocess.run(
                "module load apps/singularity; singularity exec apps/samtools.sif samtools view -bS -f 4 " + out_path +
                "_aligned.sam" + " > " + out_path + "_unmapped.bam", shell=True)
            # convert bam to fastq file
            print("Convert bam to fastq file for Import non-human reads as BAM file for " + out_path + "_unmapped.bam")
            subprocess.run(
                "module load apps/singularity; singularity exec apps/samtools.sif samtools bam2fq " + out_path +
                "_unmapped.bam > " + out_path + "_unmapped.fastq", shell=True)
            # remove intermediary file
            os.remove(out_path + "_unmapped.bam")
    print("--------------------------")


def de_novo_assembly(data_dir, out_dir, run_id):
    # --meta enables mode for metagenome/uneven coverage assembly, designed for highly non-uniform coverage and
    # is sensitive to underrepresented sequence at low coverage (as low as 2x)
    # alternative haplotypes are collapsed
    # minimum overlap length is chosen automatically based on read length distribution
    # flye performs one polishing iteration by default
    print("--------------------------\nDE NOVO ASSEMBLY\n--------------------------")
    for data_input in glob.glob(data_dir + "/*.fq"):
        barcode = get_identifier(file=data_input, string="barcode")
        for file in glob.glob(out_dir + "/human_read_removal/" + run_id + "*" + barcode + "*.fastq"):
            if os.path.exists(out_dir + "/de_novo_assembly/" + run_id + "_" + barcode):
                print("De novo assembly already complete for " + barcode)
            else:
                print("Performing de novo assembly for " + file)
                create_directory(parent_directory=out_dir + "/de_novo_assembly", directory=run_id + "_" + barcode)
                print("Conducting de novo assembly for " + file)
                subprocess.run("module load apps/singularity; singularity exec apps/flye.sif flye --nano-raw " +
                                file + " --out-dir " + out_dir + "/de_novo_assembly/" + run_id + "_" + barcode +
                                " --meta", shell=True)
    print("--------------------------")


def mapping_assembly():
    # start with minimap2
    pass


def consensus_generation():
    # bcftools consensus generation
    pass


def variant_calling():
    # SNP-sites to identify SNPs between samples
    pass


def genetic_distance():
    # SNP-dists to calculate SNP distances
    pass


def mlst():
    # multi-locus sequence typing using MLST
    pass


def snp_based_typing():
    pass


def report_generation():
    pass


def multiqc(out_dir, run_id):
    """
    Create MultiQC report, pulling in outputs from other tools
    """
    print("--------------------------\nMULTIQC\n--------------------------")
    if os.path.exists(out_dir + "/multiqc/" + run_id):
        print("MultiQC report already generated for " + run_id)
    else:
        print("Creating MultiQC report for the analysis")
        subprocess.run(
            "module load apps/singularity; singularity exec apps/multiqc.sif python -m multiqc " + out_dir + " --outdir " +
            out_dir + "/multiqc/" + run_id, shell=True)
    print("--------------------------")


def main():
    # Install containers
    install_containers.install_tools()
    # Load singularity
    for run in os.listdir("data/run_folders"):
        if os.path.isdir("data/run_folders/" + run):
            data_dir = "data/run_folders/" + run
            for summary_file in glob.glob(data_dir + "/final_summary_*.txt"):
                run_id = get_identifier(file=summary_file, string="sample_id=").rstrip("\n")
                out_dir = "output/" + run
                # create output directory per run, and subdirectories for outputs from each tool
                create_directory(parent_directory="output", directory=run_id)
                sub_directories = ["fastqc", "pycoqc", "human_read_removal", "de_novo_assembly"]
                for i in sub_directories:
                    create_directory(parent_directory=out_dir, directory=i)
                # Conduct QC analysis
                fastqc(data_dir=data_dir, out_dir=out_dir, run_id=run_id)
                pycoqc(out_dir=out_dir, data_dir=data_dir)
                # Conduct alignment and assembly
                reference_genome = "data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna"
                human_read_removal(data_dir=data_dir, out_dir=out_dir, run_id=run_id, ref=reference_genome)
                de_novo_assembly(data_dir=data_dir, out_dir=out_dir, run_id=run_id)
                multiqc(out_dir=out_dir, run_id=run_id)


if __name__ == '__main__':
    main()
