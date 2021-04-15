import glob
import os
import re
import subprocess
import sys
import shutil
import install_containers
from install_containers import app_dictionary


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
    if os.path.exists("{}/{}".format(parent_directory, directory)):
        print("{} directory already exists".format(directory))
    else:
        os.mkdir(os.path.join(parent_directory, directory))
        print("output folder {} created".format(directory))


def fastqc(data_dir, out_dir, run_id):
    """
    FastQC analysis per run folder
    """
    print("--------------------------\nFASTQC\n--------------------------")
    for file in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=file, string="barcode")
        output_filename = "{}_{}_fastqc.html".format(run_id, barcode)
        # if output file already present, do not re-analyse
        if os.path.isfile("{}/fastqc/{}".format(out_dir, output_filename)):
            print("FastQC output for run {} already exists".format(run_id))
        else:
            print("Creating FastQC file for run {}".format(run_id))
            fastqc_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} fastqc {} -o " \
                             "{}/fastqc".format(app_dictionary["fastqc"], file, out_dir)
            subprocess.run(fastqc_command, shell=True)
    print("--------------------------")


def pycoqc(data_dir, out_dir):
    """
    Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes() function.
    Then creates a pycoQC json report per barcode.
    """
    print("--------------------------\nPYCOQC\n--------------------------")
    # If pycoQC not yet run, split the summary sequencing files according to barcodes and run pycoQC
    if not os.listdir("{}/pycoqc".format(out_dir)):
        print("pycoQC not yet run")
        split_barcodes(data_dir=data_dir, out_dir=out_dir)
        for file in glob.glob("{}/pycoqc/sequencing_summary_*".format(out_dir)):
            # get barcode name from barcode_split output file names
            barcode = file.split(".", 1)[0].split("summary_", 1)[1]
            print("Creating PycoQC json report for {}".format(file))
            pycoqc_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                             "{}/pycoqc/{}_pycoQC_output.json".format(app_dictionary["pycoqc"], file, out_dir, barcode)
            subprocess.run(pycoqc_command, shell=True)
    else:
        print("Directory not empty - pycoQC already run")
    print("--------------------------")


def split_barcodes(data_dir, out_dir):
    for file in glob.glob("{}/sequencing_summary_*".format(data_dir)):
        print("Splitting summary sequencing file {} according to barcodes".format(file))
        split_barcode_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
                                "--output_unclassified --min_barcode_percent 0.0 --summary_file {} --output_dir " \
                                "{}/pycoqc".format(app_dictionary["pycoqc"], file, out_dir)
        subprocess.run(split_barcode_command, shell=True)


def human_read_removal(data_dir, out_dir, run_id, ref):
    """
    Removes human reads from the samples by alignment to the human reference genome.
    """
    print("--------------------------\nHUMAN READ REMOVAL\n--------------------------")
    for file in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=file, string="barcode")
        out_path = "{}/human_read_removal/{}_{}".format(out_dir, run_id, barcode)
        if os.path.isfile("{}_aligned.sam".format(out_path)):
            print("{} already aligned".format(barcode))
        else:
            # align reads to human reference genome using ont-specific parameters
            print("Aligning reads to human reference genome for {}".format(file))
            minimap2_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} minimap2 -ax map-ont {} {} -o " \
                               "{}_aligned.sam".format(app_dictionary["minimap2"], ref, file, out_path)
            subprocess.run(minimap2_command, shell=True)
        if os.path.isfile("{}_unmapped.fastq".format(out_path)):
            print("Human read removal already complete for {}".format(barcode))
        else:
            # import unassigned reads from sam file and convert to fastq file
            print("Import non-human reads as fastq file for {}_aligned.sam".format(out_path))
            samtools_fastq_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} samtools fastq -f 4 " \
                                     "{}_aligned.sam -0 {}_unmapped.fastq".format(app_dictionary["samtools"], out_path,
                                                                                  out_path)
            subprocess.run(samtools_fastq_command, shell=True)
        if os.path.isfile("{}_samtools_stats.txt".format(out_path)):
            print("Samtools stats already conducted for {}".barcode)
            # remove intermediary file
            # os.remove("{}_unmapped.bam".format(out_path))
            # calculate % aligned reads to human reference genome
            samtools_stats_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                                     "samtools stats {}_aligned.sam | grep ^SN | cut -f 2- > " \
                                     "{}_samtools_stats.txt".format(app_dictionary["samtools"], out_path, out_path)
            subprocess.run(samtools_stats_command, shell=True)
    print("--------------------------")


def de_novo_assembly(data_dir, out_dir, run_id):
    # --meta enables mode for metagenome/uneven coverage assembly, designed for highly non-uniform coverage and
    # is sensitive to underrepresented sequence at low coverage (as low as 2x)
    # alternative haplotypes are collapsed
    # minimum overlap length is chosen automatically based on read length distribution
    # flye performs one polishing iteration by default
    print("--------------------------\nDE NOVO ASSEMBLY\n--------------------------")
    for data_input in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=data_input, string="barcode")
        for file in glob.glob("{}/human_read_removal/{}*{}*.fastq".format(out_dir, run_id, barcode)):
            if os.path.exists("{}/de_novo_assembly/{}_{}".format(out_dir, run_id, barcode)):
                print("De novo assembly already complete for {}".format(barcode))
            else:
                print("Performing de novo assembly for {}".format(file))
                parent_directory = "{}/de_novo_assembly".format(out_dir)
                directory = "{}_{}".format(run_id, barcode)
                create_directory(parent_directory=parent_directory, directory=directory)
                create_directory(parent_directory="{}/{}".format(parent_directory, directory),
                                 directory="pyfasta_split_contigs")
                print("Conducting de novo assembly for {}".format(file))
                flye_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} flye --nano-raw {} --out-dir " \
                               "{}/de_novo_assembly/{}_{} --meta".format(app_dictionary["flye"], file,
                                                                         out_dir, run_id, barcode)
                subprocess.run(flye_command, shell=True)
    print("--------------------------")


def mapping_assembly():
    # start with minimap2
    pass


def consensus_generation():
    # bcftools consensus generation
    pass


def split_files(data_dir, out_dir, run_id, cwd):
    print("--------------------------\nSPLIT ASSEMBLY INTO CONTIGS\n--------------------------")
    for data_input in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=data_input, string="barcode")
        assembly_directory = "{}/de_novo_assembly/{}_{}".format(out_dir, run_id, barcode)
        if not os.listdir("{}/pyfasta_split_contigs".format(assembly_directory)):
            print("Splitting assembly for {} into a file per contig".format(barcode))
            os.chdir("{}/pyfasta_split_contigs".format(assembly_directory))
            faidx_command = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} ; pwd ; pyfasta split --header " \
                            "'%(seqid)s.fasta' {}/{}/assembly.fasta".format(app_dictionary["pyfasta"], cwd,
                                                                            assembly_directory)
            subprocess.run(faidx_command, shell=True)
            os.chdir(cwd)
            for file in glob.glob("{}/pyfasta_split_contigs".format(assembly_directory)):
                append_id(filename=file, uid="{}_{}".format(run_id, barcode))
        else:
            print("Assembly already split for {}".format(barcode))
    print("--------------------------")


def append_id(filename, uid):
    name, ext = os.path.splitext(filename)
    return "{uid}_{name}{ext}".format(name=name, uid=uid, ext=ext)


def mlst(data_dir, out_dir, run_id, cwd):
    """
    Multi-locus sequence typing using MLST
    """
    # FOR THIS TO WORK I THINK I WILL HAVE TO SPLIT THE ASSEMBLY FILE INTO INDIVIDUAL CONTIGS
    print("--------------------------\nMULTI LOCUS SEQUENCE TYPING\n--------------------------")
    for data_input in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=data_input, string="barcode")
        fasta_input = "{}/de_novo_assembly/{}_{}/pyfasta_split_contigs/*.fasta".format(out_dir, run_id, barcode)
        csv_output = "{}/mlst/{}_{}.csv".format(out_dir, run_id, barcode)
        if os.path.exists(csv_output):
            print("MLST already complete for {}".format(barcode))
        else:
            print("Conducting MLST for {}".format(barcode))
            mlst_command = "sudo docker run --rm -v `pwd`:`pwd` -w `{} mlst --debug {} >> " \
                           "{}".format(app_dictionary["mlst"], fasta_input, csv_output)
            subprocess.run(mlst_command, shell=True)
    print("--------------------------")
    pass


def alignment():
    """
    Multi-fasta alignment of all contigs within
    """


def variant_calling(out_dir, run_id):
    print("--------------------------\nVARIANT CALLING\n--------------------------")
    print(out_dir)
    print(run_id)
    # SNP-sites to identify SNPs between samples
    variant_calling_command = "sudo docker run --rm -v `pwd`:`pwd` -w `{} snp-sites -m -o {} " \
                              "{}".format(app_dictionary["snp-sites"], OUTPUT_FILENAME, INPUT_FILENAME)
    subprocess.run(variant_calling_command, shell=True)
    print("--------------------------")
    pass


def genetic_distance(out_dir, run_id):
    print("--------------------------\nGENETIC DISTANCE CALCULATION\n--------------------------")
    # SNP-dists to calculate SNP distances
    genetic_distance_command = "sudo docker run --rm -v `pwd`:`pwd` -w `{} snp-dists {} > " \
                               "{}".format(app_dictionary["snp-dists"], INPUT_FILENAME, OUTPUT_FILENAME.tsv)
    subprocess.run(genetic_distance_command, shell=True)
    print("--------------------------")
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
    if os.path.exists("{}/multiqc/{}".format(out_dir, run_id)):
        print("MultiQC report already generated for {}".format(run_id))
    else:
        print("Creating MultiQC report for the analysis")
        multiqc_command = "sudo docker run --rm -v `pwd`:`pwd` -w `{} multiqc {} --outdir " \
                          "{}/multiqc/{}".format(app_dictionary["multiqc"], out_dir, out_dir, run_id)
        subprocess.run(multiqc_command, shell=True)
    print("--------------------------")


def main():
    # Install containers
    for key in app_dictionary:
        install_containers.install_tools(key)
    for run in os.listdir("data/run_folders"):
        if os.path.isdir("data/run_folders/{}".format(run)):
            data_dir = "data/run_folders/{}".format(run)
            print("--------------------------\nANALYSING RUN {}\n--------------------------".format(run))
            for summary_file in glob.glob("{}/final_summary_*.txt".format(data_dir)):
                run_id = get_identifier(file=summary_file, string="sample_id=").rstrip("\n")
                out_dir = "output/{}".format(run)
                # create output directory per run, and subdirectories for outputs from each tool
                create_directory(parent_directory="output", directory=run_id)
                sub_directories = ["fastqc", "pycoqc", "human_read_removal", "de_novo_assembly", "mlst", "snp-sites",
                                   "snp-dists"]
                for i in sub_directories:
                    create_directory(parent_directory=out_dir, directory=i)
                # Conduct QC analysis
                fastqc(data_dir=data_dir, out_dir=out_dir, run_id=run_id)
                pycoqc(out_dir=out_dir, data_dir=data_dir)
                # Conduct alignment and assembly
                reference_genome = "data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna"
                human_read_removal(data_dir=data_dir, out_dir=out_dir, run_id=run_id, ref=reference_genome)
                de_novo_assembly(data_dir=data_dir, out_dir=out_dir, run_id=run_id)
                split_files(data_dir=data_dir, out_dir=out_dir, run_id=run_id, cwd=os.getcwd())
                # mlst(data_dir=data_dir, out_dir=out_dir, run_id=run_id, cwd=cwd)
                # variant_calling(out_dir=out_dir, run_id=run_id)
                # genetic_distance(out_dir=out_dir, run_id=run_id)
                # multiqc(out_dir=out_dir, run_id=run_id)


if __name__ == '__main__':
    main()
