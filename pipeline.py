import glob
import os
import subprocess
import install_containers


def get_identifier(runfolder, file, string):
    """
    Gets the run id from the run's final_summary file
    """
    with open("data/run_folders/" + runfolder + "/" + file, "rt") as file_contents:
        lines = file_contents.readlines()
        for line in lines:
            if line.__contains__(string):
                identifier = (line.split("=", 1))[0]
    return identifier


def create_directory(parent_directory, directory_list):
    """
    Create output directories
    """
    for i in directory_list:
        if os.path.exists(parent_directory + i):
            print(i + "directory already exists")
        else:
            os.mkdir(os.path.join(parent_directory, i))
            print("output folder " + i + "created")


def fastqc(run_folder, run_id):
    """
    FastQC analysis per run folder
    """
    output_filename = run_id + "_24hrs_fastqc.html"
    # if output file already present, do not re-analyse
    if os.path.isfile("output/" + run_folder + "/fastqc/" + output_filename):
        print("FastQC output for run " + run_id + " already exists")
    else:
        print("Creating fastqc file for run " + run_id)
        subprocess.run("singularity exec apps/fastqc.sif fastqc data/run_folders/" + run_folder +
                       "_24hrs.fq -o output/fastqc", shell=True)


def pycoqc(run_folder):
    """
    Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes() function.
    Then creates a pycoQC json report per barcode.
    """
    # If barcodes not yet split, split the summary sequencing files according to barcodes
    if not os.listdir("output/" + run_folder + "/pycoqc"):
        print("Summary sequencing files not yet split")
        split_barcodes(run_folder)
    else:
        print("Directory not empty - barcodes already split")

    for file in glob.glob("output/" + run_folder + "/pycoqc/sequencing_summary_*"):
        # get barcode name from barcode_split output file names
        barcode = subprocess.run("echo " + file + " | cut -d '_' -f 3 | cut -d '.' -f 1")
        print("Creating pycoQC json report for " + file)
        subprocess.run(
            "singularity exec apps/pycoqc.sif pycoQC -f " + file +
            " --json_outfile output/" + run_folder + "/pycoqc/" + barcode + "_pycoQC_output.json", shell=True)


def split_barcodes(run_folder):
    for file in glob.glob("data/run_folders/" + run_folder + "/sequencing_summary_*"):
        print("Splitting summary sequencing file " + file + " according to barcodes")
        subprocess.run(
            "singularity exec apps/pycoqc.sif Barcode_split --output_unclassified --min_barcode_percent 0.0 " +
            "--summary_file " + "data/run_folders/" + run_folder + "/" + file + " --output_dir output/pycoqc/",
            shell=True)


def human_read_removal(ref, run_folder, run_id):
    """
    Removes human reads from the samples by alignment to the human reference genome.

    """
    for file in glob.glob("data/run_folders/" + run_folder + "/*.fq"):
        barcode = get_identifier(runfolder=run_folder, file=file, string="barcode")
        # align reads to human reference genome
        print("Aligning reads to human reference genome for " + file)
        subprocess.run(
            "singularity exec apps/minimap2.sif minimap2 -ax map-ont " + ref + file +
            " > output/human_read_removal/" + run_id + "_" + barcode + "_aligned.sam", shell=True)
        # export unassigned reads to bam file with samtools
        # singularity shell apps/samtools -c "samtools view -f 4 file.bam > unmapped.sam"

        # align reads using ont-specific parameters

        # minimap2 for mapping alignment, bcftools consensus generation, SNP-sites to identify SNPs between samples
        # multi-locus sequence typing using MLST, and SNP-dists to calculate SNP distances.1


def multiqc():
    """
    Create MultiQC report, pulling in outputs from other tools
    """
    print("Creating MultiQC report for the analysis")
    subprocess.run("singularity exec apps/multiqc.sif python -m multiqc output --outdir output/multiqc", shell=True)


def main():
    # Install containers
    install_containers.install_tools()
    # Load singularity
    subprocess.run("module load apps/singularity", shell=True)
    for directory in os.listdir("data/run_folders"):
        if os.path.isdir("data/run_folders/" + directory):
            run_id = get_identifier(runfolder=dir, file="final_summary_*.txt", string="sample_id=")
            # create output directory per run, and subdirectories for outputs from each tool
            create_directory(parent_directory="output", directory_list=run_id)
            sub_directories = ["fastqc", "pycoqc", "human_read_removal"]
            create_directory(parent_directory="output/" + run_id, directory_list=sub_directories)
            # Conduct fastQC analysis
            fastqc(run_folder=directory, run_id=run_id)
            # Conduct pycoQC analysis
            pycoqc(run_folder=directory)
            reference_genome = "data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna"
            human_read_removal(ref=reference_genome, run_folder=dir, run_id=run_id)


if __name__ == '__main__':
    main()
