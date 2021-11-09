"""
Pipeline
"""
import glob
import os
import re
import subprocess
import sys
import shutil
import install_containers
import time
from install_containers import app_dictionary
import get_reference_sequences


def get_identifier(file, string):
    """
    Gets the run id from the run's final_summary file
    """
    if ".fq" in file:
        with open(file, "rt") as file_contents:
            line = file_contents.readline()
            identifier = (line.split("=", 1))[1].split()[-1].split("=")[-1]
    else:
        with open(file, "rt") as file_contents:
            lines = file_contents.readlines()
            for line in lines:
                if line.__contains__(string):
                    identifier = (line.split("=", 1))[1].split()[-1].split("=")[-1]
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
            print("FastQC output for run {} {} already exists".format(run_id, barcode))
        else:
            print("Creating FastQC file for run {}".format(run_id))
            fastqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} fastqc {} -o " \
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
            pycoqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                             "{}/pycoqc/{}_pycoQC_output.json".format(app_dictionary["pycoqc"], file, out_dir, barcode)
            subprocess.run(pycoqc_command, shell=True)
    else:
        print("Directory not empty - pycoQC already run")
    print("--------------------------")


def split_barcodes(data_dir, out_dir):
    """
    Split barcodes
    """
    for file in glob.glob("{}/sequencing_summary_*".format(data_dir)):
        print("Splitting summary sequencing file {} according to barcodes".format(file))
        split_barcode_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
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
            minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} minimap2 -ax map-ont {} {} -o " \
                               "{}_aligned.sam".format(app_dictionary["minimap2"], ref, file, out_path)
            subprocess.run(minimap2_command, shell=True)
        if os.path.isfile("{}_unmapped.fastq".format(out_path)):
            print("Human read removal already complete for {}".format(barcode))
        else:
            # import unassigned reads from sam file and convert to fastq file
            print("Import non-human reads as fastq file for {}_aligned.sam".format(out_path))
            samtools_fastq_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} samtools fastq -f 4 " \
                                     "{}_aligned.sam -0 {}_unmapped.fastq".format(app_dictionary["samtools"], out_path,
                                                                                  out_path)
            subprocess.run(samtools_fastq_command, shell=True)
        if os.path.isfile("{}_samtools_stats.txt".format(out_path)):
            print("Samtools stats already conducted for {}".format(barcode))
            # remove intermediary file
            # os.remove("{}_unmapped.bam".format(out_path))
            # calculate % aligned reads to human reference genome
        else:
            print("Run samtools stats for {}".format(barcode))
            samtools_stats_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                                     "samtools stats {}_aligned.sam | grep ^SN | cut -f 2- > " \
                                     "{}_samtools_stats.txt".format(app_dictionary["samtools"], out_path, out_path)
            subprocess.run(samtools_stats_command, shell=True)
    print("--------------------------")

def centrifuge_index(cwd):
    """
    Download centrifuge database and create index
    """
    data_dir = "data"
    directory = "centrifuge_index"
    create_directory(data_dir, directory)
    centrifuge_dir = "{}/{}/{}".format(cwd, data_dir, directory)
    os.chdir(centrifuge_dir)

    print("--------------------------\nDOWNLOAD CENTRIFUGE DATABASE\n--------------------------")
    filelist = ['names.dmp', 'nodes.dmp']
    if all([os.path.isfile("{}/{}/{}/taxonomy/{}".format(cwd, data_dir, directory, file)) for file in filelist]):
        print("Centrifuge database already downloaded")
    else:
        # masks low complexity regions using -m
        # only downloads species from the bacteria domain, matching the listed taxa IDs
        centrifuge_taxonomy = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} " \
                              "centrifuge-download -o taxonomy taxonomy".format(app_dictionary["centrifuge"])
        centrifuge_library = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} centrifuge-download -o library -m " \
                             "-d 'bacteria' -t 1352,354276,562,573,287,615,1280,1352 refseq > " \
                             "seqid2taxid.map".format(app_dictionary["centrifuge"])

        print("Downloading taxonomy...")
        subprocess.run(centrifuge_taxonomy, shell=True)
        print("Downloading reference sequences...")
        subprocess.run(centrifuge_library, shell=True)


    print("--------------------------\nCREATE CENTRIFUGE INDEX\n--------------------------")
    if os.path.isfile("{}/input-sequences.fna".format(centrifuge_dir)):
        print("Index file already created")
    else:
        print("Concatenating input sequences...")
        index = "cat library/*/*.fna > input-sequences.fna"
        subprocess.run(index, shell=True)

    print("Building index...")
    centrifuge_build = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} centrifuge-build -p 4 --conversion-table " \
                       "seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp " \
                       "input-sequences.fna index".format(app_dictionary["centrifuge"])
    subprocess.run(centrifuge_build, shell = True)
    # need to add some way to parse the centrifuge output and use this to filter the input reads to return only those
    # that returned a classification from centrifuge (i.e. are from the species of interest)


def classification(cwd):
    """
    Centrifuge classification
    Call centrifuge_index to create the index if not already created
    Call centrifuge
    """
    # create centrifuge index
    centrifuge_index(cwd)
    # run centrifuge-inspect to output a summary of the index used for classification
    centrifuge_inspect = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} centrifuge-inspect -s " \
                         "index".format(app_dictionary["centrifuge"])
    os.chdir(cwd)
    subprocess.run(centrifuge_inspect, shell = True)



def de_novo_assembly(data_dir, out_dir, run_id, cwd):
    """
    Runs using flye.
    """
    # --meta enables mode for metagenome/uneven coverage assembly, designed for highly non-uniform coverage and
    # is sensitive to underrepresented sequence at low coverage (as low as 2x)
    # alternative haplotypes are collapsed
    # minimum overlap length is chosen automatically based on read length distribution
    # flye performs one polishing iteration by default
    print("--------------------------\nDE NOVO ASSEMBLY\n--------------------------")
    for data_input in glob.glob("{}/{}/*.fq".format(cwd, data_dir)):
        barcode = get_identifier(file=data_input, string="barcode")
        for file in glob.glob("{}/human_read_removal/{}*{}*.fastq".format(out_dir, run_id, barcode)):
            assembly_output = "assembly_{}_{}".format(run_id, barcode)
            parent_directory = "{}/de_novo_assembly".format(out_dir)
            directory = "{}_{}".format(run_id, barcode)
            if os.path.isfile("{}/{}/{}.fasta".format(parent_directory, directory, assembly_output)):
                print("De novo assembly already complete for {}".format(barcode))
            else:
                print("Performing FLYE de novo assembly for {}".format(file))
                create_directory(parent_directory=parent_directory, directory=directory)
                #create_directory(parent_directory="{}/{}".format(parent_directory, directory),
                #                 directory="pyfasta_split_contigs")
                print("Conducting de novo assembly for {}".format(file))
                assembly_outdir = parent_directory + "/" + directory
                flye_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} flye --nano-raw {} --out-dir {} " \
                               "--meta".format(app_dictionary["flye"], file, assembly_outdir)
                subprocess.run(flye_command, shell=True)
                # ADD LINE TO RENAME OUTPUT TO ADD IN THE RUN ID AND BARCODE
            assembly_proc_out(parent_directory, assembly_output, directory, barcode)
    print("--------------------------")


def assembly_proc_out(parent_directory, assembly_output, directory, barcode):
    """
    Save only sequences from mlst output to file.
    ADD FUNCTIONALITY TO ONLY KEEP CONTIGS THAT HAVE COVERAGE OVER
    """
    in_fasta = "{}/{}/{}.fasta".format(parent_directory, directory, assembly_output)
    searchquery = '>'
    if os.path.isfile("{}/{}/processed_*".format(parent_directory, directory)):
        print("Assembly output already processed for {}".format(barcode))
    else:
        print("Processing assembly output for {}".format(barcode))
        # ADD IN CODE TO PARSE THE assembly_info.txt FILE AND ONLY GET LINES THAT MATCH THE CONTIGS FROM THE TEXT FILE
        # WITH COVERAGE OVER 60X
        with open(in_fasta) as infile:
                lines = infile.readlines()
                for i, line in enumerate(lines):
                    if line.startswith(searchquery):
                        contig = line.split(' ')[0].split('>')[1]
                        output = "{}/{}/processed_{}_{}.fasta".format(parent_directory, directory,
                                                                      assembly_output, contig)
                        with open(output, 'a') as outfile:
                            outfile.write(line)
                            outfile.write(lines[i + 1])


def mlst(data_dir, out_dir, run_id):
    """
    Multi-locus sequence typing using MLSTcheck.
    Databases are bundled with the container so this doesn't require a separate download.
    """
    # FOR THIS TO WORK I THINK I WILL HAVE TO SPLIT THE ASSEMBLY FILE INTO INDIVIDUAL CONTIGS
    print("--------------------------\nMULTI LOCUS SEQUENCE TYPING\n--------------------------")
    for data_input in glob.glob("{}/*.fq".format(data_dir)):
        barcode = get_identifier(file=data_input, string="barcode")
        fasta_string = "{}/de_novo_assembly/{}_{}/processed_assembly_{}_{}*".format(out_dir, run_id, barcode, run_id,
                                                                                    barcode)
        csv_output = "{}/mlst/{}_{}.csv".format(out_dir, run_id, barcode)
        parent_directory = "{}/de_novo_assembly".format(out_dir)
        directory = "mlst"
        if os.path.exists(csv_output):
            print("MLST already complete for {}".format(barcode))
        else:
            print("Conducting MLST for {}".format(barcode))
            create_directory(parent_directory=parent_directory, directory=directory)
            mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} get_sequence_type {} -o " \
                           "{}/{}".format(app_dictionary["mlstcheck"], fasta_string, out_dir, directory)
            print(mlst_command)
            #subprocess.run(mlst_command, shell=True)
    print("--------------------------")
    pass


def variant_calling(out_dir, run_id):
    """
    Identify SNPs between samples
    """
    print("--------------------------\nVARIANT CALLING\n--------------------------")
    print(out_dir)
    print(run_id)
    # SNP-sites to identify SNPs between samples
    variant_calling_command = "docker run --rm -v `pwd`:`pwd` -w `{} snp-sites -m -o {} " \
                              "{}".format(app_dictionary["snp-sites"], OUTPUT_FILENAME, INPUT_FILENAME)
    subprocess.run(variant_calling_command, shell=True)
    print("--------------------------")
    pass


def genetic_distance():
    """
    Calculate SNP distances
    """
    print("--------------------------\nGENETIC DISTANCE CALCULATION\n--------------------------")
    # SNP-dists to calculate SNP distances
    genetic_distance_command = "docker run --rm -v `pwd`:`pwd` -w `{} snp-dists {} > " \
                               "{}".format(app_dictionary["snp-dists"], INPUT_FILENAME, OUTPUT_FILENAME.tsv)
    subprocess.run(genetic_distance_command, shell=True)
    print("--------------------------")
    pass


def report_generation():
    """
    Report containing information about the identified sequences.
    """
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
        multiqc_command = "docker run --rm -v `pwd`:`pwd` -w `{} multiqc {} --outdir " \
                          "{}/multiqc/{}".format(app_dictionary["multiqc"], out_dir, out_dir, run_id)
        subprocess.run(multiqc_command, shell=True)
    print("--------------------------")


def main():
    cwd = os.getcwd()

    # Get reference sequences
    out_dir = "{}/data/reference_sequences/".format(cwd)
    create_directory("data", "reference_sequences")
    get_reference_sequences.download_sequences(get_reference_sequences.refseq_dict, out_dir)

    # Install containers
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])
    for run in os.listdir("data/run_folders"):
        if os.path.isdir("data/run_folders/{}".format(run)):
            data_dir = "data/run_folders/{}".format(run)
            print("--------------------------\nANALYSING RUN {}\n--------------------------".format(run))
            for summary_file in glob.glob("{}/final_summary_*.txt".format(data_dir)):
                run_id = get_identifier(file=summary_file, string="sample_id=").rstrip("\n")
                out_dir = "{}/output/{}".format(cwd, run)
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
                classification(cwd)
                de_novo_assembly(data_dir=data_dir, out_dir=out_dir, run_id=run_id, cwd=os.getcwd())
                # split_files(data_dir=data_dir, out_dir=out_dir, run_id=run_id, cwd=os.getcwd())
                mlst(data_dir=data_dir, out_dir=out_dir, run_id=run_id)
                # genetic_distance(out_dir=out_dir, run_id=run_id)
                # variant_calling(out_dir=out_dir, run_id=run_id)
                # multiqc(out_dir=out_dir, run_id=run_id)
                # report_generation()


if __name__ == '__main__':
    main()


# def split_files(data_dir, out_dir, run_id, cwd):
#     print("--------------------------\nSPLIT ASSEMBLY INTO CONTIGS\n--------------------------")
#     for data_input in glob.glob("{}/*.fq".format(data_dir)):
#         barcode = get_identifier(file=data_input, string="barcode")
#         assembly_directory = "{}/de_novo_assembly/{}_{}".format(out_dir, run_id, barcode)
#         if (not os.listdir("{}/pyfasta_split_contigs".format(assembly_directory))) and \
#                 (os.path.isfile("{}/assembly.fasta".format(assembly_directory))):
#             print("Splitting assembly for {} into a file per contig".format(barcode))
#             os.chdir("{}/pyfasta_split_contigs".format(assembly_directory))
#             faidx_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} ; pwd ; pyfasta split --header " \
#                             "'%(seqid)s.fasta' {}/{}/assembly.fasta".format(app_dictionary["pyfasta"], cwd,
#                                                                             assembly_directory)
#             subprocess.run(faidx_command, shell=True)
#             os.chdir(cwd)
#             for file in glob.glob("{}/pyfasta_split_contigs".format(assembly_directory)):
#                 append_id(filename=file, uid="{}_{}".format(run_id, barcode))
#         else:
#             print("Assembly already split for {}".format(barcode))
#     print("--------------------------")
#
#
# def append_id(filename, uid):
#     name, ext = os.path.splitext(filename)
#     return "{uid}_{name}{ext}".format(name=name, uid=uid, ext=ext)

#
# def mapping_assembly():
#     # start with minimap2
#     pass
#
#
# def consensus_generation():
#     # bcftools consensus generation
#     pass
