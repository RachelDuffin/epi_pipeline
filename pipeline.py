import glob
import sys
import os
import subprocess
import re
import install_containers

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
        print("output folder " + directory + "created")


def fastqc(run_folder, run_id):
    """
    FastQC analysis per run folder
    """
    print("--------------------------\nFASTQC\n--------------------------")
    output_filename = run_id + "_24hrs_fastqc.html"
    # if output file already present, do not re-analyse
    if os.path.isfile("output/" + run_folder + "/fastqc/" + output_filename):
        print("FastQC output for run " + run_id + " already exists")
    else:
        print("Creating fastqc file for run " + run_id)
        for file in glob.glob("data/run_folders/" + run_folder + "*_24hrs.fq"):
            barcode = get_identifier(file=file, string="barcode")
            subprocess.run("module load apps/singularity; singularity exec apps/fastqc.sif fastqc data/run_folders/" + run_folder +
                       barcode + "_24hrs.fq -o output/" + run_folder + "/fastqc", shell=True)
    print("--------------------------")


def pycoqc(run_folder):
    """
    Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes() function.
    Then creates a pycoQC json report per barcode.
    """
    print("--------------------------\nPYCOQC\n--------------------------")
    # If pycoQC not yet run, split the summary sequencing files according to barcodes and run pycoQC
    if not os.listdir("output/" + run_folder + "/pycoqc"):
        print("pycoQC not yet run")
        split_barcodes(run_folder)
        for file in glob.glob("output/" + run_folder + "/pycoqc/sequencing_summary_*"):
            # get barcode name from barcode_split output file names
            barcode = file.split(".", 1)[0].split("summary_", 1)[1]
            print("Creating pycoQC json report for " + file)
            subprocess.run(
                "module load apps/singularity; singularity exec apps/pycoqc.sif pycoQC -f " + file +
                " --json_outfile output/" + run_folder + "/pycoqc/" + barcode + "_pycoQC_output.json", shell=True)
    else:
        print("Directory not empty - pycoQC already run")
    print("--------------------------")


def split_barcodes(run_folder):
    for file in glob.glob("data/run_folders/" + run_folder + "/sequencing_summary_*"):
        print("Splitting summary sequencing file " + file + " according to barcodes")
        subprocess.run(
            "module load apps/singularity; singularity exec apps/pycoqc.sif Barcode_split --output_unclassified " +
            "--min_barcode_percent 0.0 --summary_file " + file + " --output_dir output/" + run_folder + "/pycoqc",
            shell=True)


def human_read_removal(ref, run_folder, run_id):
    """
    Removes human reads from the samples by alignment to the human reference genome.

    """
    print("--------------------------\nHUMAN READ REMOVAL\n--------------------------")
    for file in glob.glob("data/run_folders/" + run_folder + "/*.fq"):
        barcode = get_identifier(file=file, string="barcode")
        if os.path.isfile("output/" + run_folder + "/human_read_removal/" + run_folder + "_" + barcode + "_aligned.sam"):
            print(barcode + " already aligned")
        else:
            # align reads to human reference genome using ont-specific parameters
            print("Aligning reads to human reference genome for " + file)
            subprocess.run(
                "module load apps/singularity; singularity exec apps/minimap2.sif minimap2 -ax map-ont " + ref + " " +
                file + " > output/" + run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_aligned.sam",
                shell=True)
        # import SAM to BAM as @SQ lines present in header, only import unassigned reads (non-human)
        print("Import non-human reads as BAM file for output/" + run_folder + "/human_read_removal/" + run_id + "_"
              + barcode + "_aligned.sam")
        subprocess.run(
            "module load apps/singularity; singularity exec apps/samtools.sif samtools view -bS -f 4 " + "output/" +
            run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_aligned.sam" + " > " + "output/" +
            run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_unmapped.bam", shell=True)
        # convert bam to fastq file
        print("Convert bam to fastq file for Import non-human reads as BAM file for output/" +
            run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_unmapped.bam")
        subprocess.run(
            "module load apps/singularity; singularity exec apps/samtools.sif samtools bam2fq output/" +
            run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_unmapped.bam > output/" +
            run_folder + "/human_read_removal/" + run_id + "_" + barcode + "_unmapped.fastq", shell=True)
    print("--------------------------")

        # minimap2 for mapping alignment, bcftools consensus generation, SNP-sites to identify SNPs between samples
        # multi-locus sequence typing using MLST, and SNP-dists to calculate SNP distances.1

def de_novo_assembly():
    pass

def mapping_assembly():
    pass

def consensus_generation():
    pass

def variant_calling():
    pass

def genetic_distance():
    pass

def mlst():
    pass

def snp_based_typing():
    pass

def report_generation():
    pass

def multiqc():
    """
    Create MultiQC report, pulling in outputs from other tools
    """
    print("--------------------------\nMULTIQC\n--------------------------")
    print("Creating MultiQC report for the analysis")
    subprocess.run("module load apps/singularity; singularity exec apps/multiqc.sif python -m multiqc output --outdir output/multiqc", shell=True)
    print("--------------------------")

def main():
    # Install containers
    install_containers.install_tools()
    # Load singularity
    for run in os.listdir("data/run_folders"):
        if os.path.isdir("data/run_folders/" + run):
            for summary_file in glob.glob("data/run_folders/" + run + "/final_summary_*.txt"):
                run_id = get_identifier(file=summary_file, string="sample_id=").rstrip("\n")
                # create output directory per run, and subdirectories for outputs from each tool
                create_directory(parent_directory="output", directory=run_id)
                sub_directories = ["fastqc", "pycoqc", "human_read_removal"]
                for i in sub_directories:
                    create_directory(parent_directory="output/" + run_id, directory=i)
                # Conduct fastQC analysis
                fastqc(run_folder=run, run_id=run_id)
                # Conduct pycoQC analysis
                pycoqc(run_folder=run)
                reference_genome = "data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna"
                human_read_removal(ref=reference_genome, run_folder=run, run_id=run_id)


if __name__ == '__main__':
    main()
