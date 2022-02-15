"""
Pipeline
"""
import glob
import os
import re
from subprocess import run, Popen, PIPE, STDOUT
import shlex
import sys
import shutil
import install_containers
import time
from config import app_dictionary
import get_reference_sequences
import logging
import config
import argparse

class Pipeline:
    """Loop through and process runfolders in a given directory.
    A single class instance is required to process all runfolders.

    Methods defined here:
    loop_through_runs()

    """

    def __init__(self, run, runfolder_dir, reference_dir, index_dir, out_dir, ref_gen):
        """
        The constructor for Pipeline class
        """
        self.now = ""  # Stores time stamp for class instance, used in log file.

        # specify data inputs paths
        self.reference_genome = ref_gen
        self.app_dictionary = app_dictionary
        self.run = run
        self.runfolder_dir = runfolder_dir
        self.runfolder = "{}/{}".format(self.runfolder_dir, self.run)
        self.sequencing_summary_file = glob.glob("{}/final_summary_*.txt".format(self.runfolder))
        self.run_id = self.get_identifier(self.sequencing_summary_file, "sample_id=").rstrip("\n")
        self.fastq_path = "{}/*.fq".format(self.runfolder)
        self.reference_dir = reference_dir
        self.index_dir = index_dir

        # specify outputs
        self.out_dir = out_dir
        self.runfolder_out = "{}/{}".format(out_dir, self.run)
        self.centrifuge_index_summary = "{}/centrifuge_index_summary.txt".format(self.runfolder_out)
        self.out_subdirs = ["QC", "human_read_removal", "de_novo_assembly", "mlst", "snp-sites",
                           "snp-dists", "centrifuge", "split_assemblies", "multiqc"]

        self.hr_rem_dir = "{}/human_read_removal".format(self.runfolder_out)
        self.class_dir = "{}/centrifuge".format(self.runfolder_out)

        self.index = "{}//p_compressed+h+v".format(self.index_dir)
        self.assembly_dir = "{}/de_novo_assembly".format(self.runfolder_out)
        self.assembly_split_dir = "{}/split_assemblies".format(self.runfolder_out)
        self.mlst_dir = "{}/mlst".format(self.runfolder_out)
        self.qc_dir = "{}/QC".format(self.runfolder_out)
        self.multiqc_dir = "{}/multiqc".format(self.runfolder_out)

        self.logfile = "{}/{}_logfile.txt".format(self.runfolder_out, self.run_id)


    def set_off_analysis(self):
        """
        Per-run analysis
        """

        if os.path.isdir(self.runfolder):
            # create run output dir if it doesn't already exist
            self.create_directory(self.runfolder_out)

            # create logfile if it doesn't already exist
            if not os.path.isfile(self.logfile):
                open(self.logfile, 'w')

            # create output directory for run
            self.create_directory(self.out_dir)

            # create subdirectories for outputs from each tool
            for subdir in self.out_subdirs:
                out_subdir = "{}/{}".format(self.runfolder_out, subdir)
                self.create_directory(out_subdir)

            # get any reference material that isn't already present
            self.get_centrifuge_index()
            get_reference_sequences.download_sequences(config.refseq_dict, self.reference_dir)

            self.logger().info("ANALYSING RUN {}".format(self.run))
            # run pycoqc for whole run
            self.pycoqc()

            for fastq in glob.glob(self.fastq_path):

                barcode = self.get_identifier(file_list=[fastq], string="barcode")
                run_barcode = "{}_{}".format(self.run_id, barcode)

                self.fastqc(fastq, run_barcode)
                unmapped_fastq = self.human_read_removal(fastq, run_barcode)
                self.classification(run_barcode, unmapped_fastq)
                assembly, assembly_outdir = self.de_novo_assembly(run_barcode, unmapped_fastq)
                # get dictionary of contigs that pass coverage requirements
                # passing_contigs = self.get_passing_contigs(run_barcode, assembly_outdir)

                #assembly_dict = self.assembly_split(run_barcode, assembly, assembly_outdir, passing_contigs)

                self.mlst(run_barcode, assembly)

                    # genetic_distance(out_dir=out_dir, run_id=run_id)
                    # variant_calling(out_dir=out_dir, run_id=run_id)
            self.multiqc(run_barcode)
                    # report_generation()

    def get_passing_contigs(self, run_barcode, assembly_outdir):
        """
        Parse assembly info file to get a list of contigs that pass minimum coverage requirements (60X)
        """
        contig_list = []
        assembly_info = "{}/assembly_info.txt".format(assembly_outdir)

        if os.path.isfile(assembly_info):
            self.logger().info("PARSE ASSEMBLY INFO: Parsing assembly info file for {}".format(run_barcode))
            with open(assembly_info, "r") as file:
                for line in file.readlines():
                    if not line.startswith("#seq_name"):
                        line = line.split("\t")
                        coverage = int(line[2])
                        if coverage >= 60:
                            contig_list.append(line[0])
        else:
            self.logger().warning("PARSE ASSEMBLY INFO: Assembly info file does not exist for {}".format(run_barcode))

        return contig_list

    @staticmethod
    def get_identifier(file_list, string):
        """
        Gets the run id from an input file (sequencing summary file or fastq file)
        """
        for file in file_list:
            # fastq file
            if ".fq" in file:
                with open(file, "rt") as file_contents:
                    line = file_contents.readline()
                    identifier = (line.split("=", 1))[1].split()[-1].split("=")[-1]
            # sequencing summary file
            else:
                with open(file, "rt") as file_contents:
                    lines = file_contents.readlines()
                    for line in lines:
                        if line.__contains__(string):
                            identifier = (line.split("=", 1))[1].split()[-1].split("=")[-1]
        return identifier


    def pycoqc(self):
        """
        Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes() function.
        Then creates a pycoQC json report per barcode.
        """
        # If pycoQC not yet run, split the summary sequencing files according to barcodes and run pycoQC
        if not os.path.isfile("{}/{}_pycoQC_output.json".format(self.qc_dir, self.run_id)):
            # specify output directories
            self.split_barcodes()

            for file in glob.glob("{}/sequencing_summary_*".format(self.runfolder)):
                self.logger().info("PYCOQC - Creating PycoQC json report for {}. ".format(self.run_id))
                self.logger().info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
                pycoqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                                 "{}/{}_pycoQC_output.json".format(self.app_dictionary["pycoqc"], file,
                                                                          self.qc_dir, self.run_id)
                self.run_command(pycoqc_command, "PYCOQC")
        else:
            self.logger().info("PYCOQC - Directory not empty - pycoQC already run for {}".format(self.run_id))


    def split_barcodes(self):
        """
        Split barcodes
        """
        if not glob.glob("{}/sequencing_summary_*".format(self.qc_dir)):

            for file in glob.glob("{}/sequencing_summary_*".format(self.runfolder)):
                self.logger().info("PYCOQC BARCODE SPLIT: Splitting summary sequencing file {} "
                                   "according to barcodes".format(file))
                self.logger().info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
                split_barcode_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
                                        "--output_unclassified --min_barcode_percent 0.0 --summary_file {} " \
                                        "--output_dir {}".format(self.app_dictionary["pycoqc"], file, self.qc_dir)
                self.run_command(split_barcode_command, "PYCOQC BARCODE SPLIT")
        else:
            self.logger().warning("PYCOQC BARCODE SPLIT: Summary "
                                  "sequencing file already split for {} ".format(self.run_id))


    def fastqc(self, fastq, run_barcode):
        """
        FastQC analysis per barcode
        """
        fastqc_output = run_barcode + "_fastqc.html"
        # if output file already present, do not re-analyse
        if os.path.isfile("{}/{}".format(self.qc_dir, fastqc_output)):
            self.logger().warning("FASTQC: Output for {} already exists".format(run_barcode))
        else:
            self.logger().info("FASTQC: Creating FastQC file for {}".format(run_barcode))
            self.logger().info("FASTQC - version {}".format(self.app_dictionary["fastqc"]))
            fastqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                             "fastqc {} -o {}".format(self.app_dictionary["fastqc"], fastq, self.qc_dir)
            self.run_command(fastqc_command, "FASTQC")


    def human_read_removal(self, fastq, run_barcode):
        """
        Removes human reads from the samples by alignment to the human reference genome.
        """
        hr_rem_outpath = "{}/{}".format(self.hr_rem_dir, run_barcode)
        unmapped_fastq = "{}_unmapped.fastq".format(hr_rem_outpath)
        sam_output = "{}_aligned.sam".format(hr_rem_outpath)
        samtools_stats_file = "{}/{}_samtools_stats.txt".format(self.qc_dir, run_barcode)

        if os.path.isfile(sam_output):
            self.logger().warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            # align reads to human reference genome using ont-specific parameters
            self.logger().info("HR REMOVAL: Aligning reads to human reference genome for {}".format(fastq))
            self.logger().info("MINIMAP2 - version {}".format(self.app_dictionary["minimap2"]))
            minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} minimap2 -ax " \
                               "map-ont {} {} -o {}".format(self.app_dictionary["minimap2"],
                                                            self.reference_genome, fastq, sam_output)
            self.run_command(minimap2_command, "MINIMAP2")

        if os.path.isfile(samtools_stats_file):
            self.logger().warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            # convert unassigned (non-human) reads from sam file to fastq file
            self.logger().info("SAMTOOLS FASTQ: Convert non-human from sam "
                               "to fastq for {}_unmapped.sam".format(hr_rem_outpath))
            self.logger().info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))
            samtools_fastq_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} " \
                                     "samtools fastq -f 4 {} -0 {}".format(self.app_dictionary["samtools"],
                                                                            sam_output, unmapped_fastq)
            self.run_command(samtools_fastq_command, "SAMTOOLS FASTQ")

        if os.path.isfile(samtools_stats_file):
            self.logger().warning("SAMTOOLS STATS: Samtools stats already conducted for {}".format(run_barcode))
                # remove intermediary file
                # os.remove("{}_unmapped.bam".format(out_path))
                # calculate % aligned reads to human reference genome
        else:
            self.logger().info("SAMTOOLS STATS: Run samtools stats for {}".format(run_barcode))
            self.logger().info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))
            samtools_stats_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} samtools stats " \
                                     "{} | grep ^SN | cut -f 2- > {}".format(self.app_dictionary["samtools"],
                                                                             sam_output, samtools_stats_file)
            self.run_command(samtools_stats_command, "SAMTOOLS STATS")

        return unmapped_fastq


    # def concatenate_input_sequences(self):
    #     """
    #     Concatenate input sequences into one fna file
    #     """
    #     if os.path.isfile("{}/input-sequences.fna".format(self.index_dir)):
    #         print("Library already concatenated")
    #     else:
    #         print("Concatenating library...")
    #         index = "cat library/*/*.fna > input-sequences.fna"
    #         subprocess.run(index, shell=True)


    def classification(self, run_barcode, unmapped_fastq):
        """
        Centrifuge classification
        Create centrifuge index summary file, and run centriuge
        Disregard abundance - there is a bug in centrifuge that calculates the abundance as 0 for all taxIDs
        """
        centrifuge_results = "{}/{}_results.tsv".format(self.class_dir, run_barcode)
        centrifuge_summary = "{}/{}_summary.tsv".format(self.class_dir, run_barcode)
        centrifuge_kreport = "{}/{}_kreport.tsv".format(self.class_dir, run_barcode)

        # Create centrifuge index summary
        if os.path.isfile(self.centrifuge_index_summary):
            self.logger().warning("CENTRIFUGE: Centrifuge index summary file already created")
        else:
            self.logger().info("CENTRIFUGE: Creating centrifuge index summary file for {}".format(run_barcode))
            self.logger().info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            # run centrifuge-inspect to output a summary of the index used for classification
            centrifuge_inspect = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} centrifuge-inspect -s " \
                                 "{} > {}".format(self.app_dictionary["centrifuge"], self.index,
                                                  self.centrifuge_index_summary)
            self.run_command(centrifuge_inspect, "CENTRIFUGE")

        # Run centrifuge classification
        if os.path.isfile(centrifuge_results) and os.path.isfile(centrifuge_summary):
            self.logger().warning("CENTRIFUGE: Centrifuge already run for {}".format(run_barcode))
        else:
            # classify reads
            self.logger().info("CENTRIFUGE: Running centrifuge classification for {}".format(run_barcode))
            self.logger().info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            # -S is file to write classification results to, --report-file is file to write classification summary to,
            # -x is the index, -f is the input sequence, --env sets centriuge indexes environment variable so it knows
            # where to look for the index
            centrifuge_command = "docker run --env CENTRIFUGE_INDEXES={} --rm -v `pwd`:`pwd` -w `pwd` -it {} " \
                                 "centrifuge -x p_compressed+h+v -q {} -S {} " \
                                 "--report-file {}".format(self.index_dir, self.app_dictionary['centrifuge'],
                                                           unmapped_fastq, centrifuge_results, centrifuge_summary)
            self.run_command(centrifuge_command, "CENTRIFUGE")

        if os.path.isfile(centrifuge_kreport):
            self.logger().warning("CENTRIFUGE: Centrifuge-kreport already generated for {}".format(run_barcode))
        else:
            self.logger().info("CENTRIFUGE: Generating centrifuge-kreport for {}".format(run_barcode))
            self.logger().info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            kreport_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} centrifuge-kreport -x {}" \
                              " {} > {}".format(self.app_dictionary["centrifuge"], self.index,
                                                centrifuge_results, centrifuge_kreport)
            self.run_command(kreport_command, "CENTRIFUGE")

    def classification_proc_out(self):
        """
        Write function to parse fasta files and extract only those sequences that were classified to a new file
        """


    def de_novo_assembly(self, run_barcode, unmapped_fastq):
        """
        Runs using flye.
            --meta enables mode for metagenome/uneven coverage assembly, designed for highly non-uniform coverage and
            is sensitive to underrepresented sequence at low coverage (as low as 2x)
            alternative haplotypes are collapsed
            minimum overlap length is chosen automatically based on read length distribution
            flye performs one polishing iteration by default
        """
        # per-barcode output directory so intermediary files aren't overwritten in subsequent runs

        assembly_outdir = "{}/{}".format(self.assembly_dir, run_barcode)
        assembly_outfile = "{}/assembly.fasta".format(assembly_outdir)
        assembly_outfile_renamed = "{}/{}_assembly.fasta".format(assembly_outdir, run_barcode)
        assembly_log = "{}/flye.log".format(assembly_outdir)

        if os.path.isfile(assembly_log):
            self.logger().warning("FLYE: De novo assembly already complete for {}".format(run_barcode))
        else:
            self.logger().info("FLYE: Performing FLYE de novo assembly for {}".format(unmapped_fastq))
            self.logger().info("FLYE - version {}".format(self.app_dictionary["flye"]))
            self.create_directory(assembly_outdir)
            flye_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} flye --nano-raw {} --out-dir {} " \
                           "--meta".format(self.app_dictionary["flye"], unmapped_fastq, assembly_outdir)

            # rename output to add in run ID and barcode
            self.run_command(flye_command, "FLYE")

            if os.path.exists(assembly_outfile):
                try:
                    os.rename(assembly_outfile, assembly_outfile_renamed)
                except:
                    self.logger().error("FLYE: Assembly outfile {} could not be renamed".format(assembly_outfile))
            else:
                self.logger().warning("FLYE: No assembly produced.")

        return assembly_outfile_renamed, assembly_outdir


    def assembly_split(self, run_barcode, assembly, assembly_outdir, passing_contigs):
        """
        Split assembly output into a file per assembly using faidx.
        ADD FUNCTIONALITY TO ONLY KEEP CONTIGS THAT HAVE COVERAGE OVER
        """
        if os.listdir(self.assembly_split_dir):
            self.logger().warning("FAIDX: Assembly output already processed for {}".format(run_barcode))

        else:
            if os.path.exists(assembly):
                self.logger().info("FAIDX: Processing assembly output for {}".format(run_barcode))
                self.logger().info("PYFAIDX - version {}".format(self.app_dictionary["pyfaidx"]))
                faidx_command = "docker run -v `pwd`:`pwd` -v {}:{} {} /bin/sh -c 'cd {} ;" \
                                " faidx -x {}'".format(assembly_outdir, assembly_outdir, self.app_dictionary["pyfaidx"],
                                                      self.assembly_split_dir, assembly)
                self.run_command(faidx_command, "FAIDX")
            else:
                self.logger().info("FAIDX: Assembly file not present for {}".format(run_barcode))

        # create dictionary to store contig names and file paths
        assembly_dict = {}

        for file in os.listdir(self.assembly_split_dir):
            contig = file.split(".fasta")[0].rstrip("\n")
            filepath = "{}/{}".format(self.assembly_split_dir, file)
            assembly_dict[contig] = filepath

        # remove contigs not present in passing_contigs list (required to surpass coverage threshold of 60X)
        assembly_dict = self.remove_contigs(assembly_dict, passing_contigs)

        return assembly_dict


    @staticmethod
    def remove_contigs(assembly_dict, passing_contigs):
        """
        Remove contigs from assembly dictionary that are not in the list of contigs that surpass the coverage threshold
        of 60X
        """
        to_delete = []
        for item in assembly_dict:
            if item not in passing_contigs:
                to_delete.append(item)
        for item in to_delete:
            del assembly_dict[item]
        return assembly_dict


    def mlst(self, run_barcode, assembly):
        """
        Multi-locus sequence typing using MLSTcheck.
        Databases are bundled with the container so this doesn't require a separate download.
        """
        # assembly_string = ""
        # for item in assembly_dict:
        #     assembly_string = assembly_string.join(' ' + assembly_dict[item])
        # print(assembly_string)
        # #assembly_string = ' '.join(assembly_dict)

        mlst_output = "{}/{}_mlst.json".format(self.mlst_dir, run_barcode)
        #assembly_files = "{}/*".format(self.assembly_split_dir)
        if os.path.exists(mlst_output):
            self.logger().warning("MLST: mlst already complete for {}".format(run_barcode))
        else:
            self.logger().info("MLST: Conducting MLST for {}".format(run_barcode))
            self.logger().info("MLST - version {}".format(self.app_dictionary["mlst"]))

            mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` " \
                            "-it {} mlst -q --json {} {}".format(self.app_dictionary["mlst"], mlst_output,
                                                                 assembly)
            # mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` " \
            #                 "-it {} mlst -q --json {} {}".format(self.app_dictionary["mlst"], mlst_output,
            #                                                      assembly_files)
            self.run_command(mlst_command, "MLST")


    def variant_calling(self):
        """
        Identify SNPs between samples
        """
        print("--------------------------\nVARIANT CALLING\n--------------------------")
        # SNP-sites to identify SNPs between samples
        variant_calling_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} snp-sites -m -o {} " \
                                  "{}".format(self.app_dictionary["snp-sites"], OUTPUT_FILENAME, INPUT_FILENAME)
        run(variant_calling_command, shell=True)
        print("--------------------------")


    def genetic_distance(self):
        """
        Calculate SNP distances
        """
        print("--------------------------\nGENETIC DISTANCE CALCULATION\n--------------------------")
        # SNP-dists to calculate SNP distances
        genetic_distance_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} snp-dists {} > " \
                                   "{}".format(self.app_dictionary["snp-dists"], INPUT_FILENAME, OUTPUT_FILENAME.tsv)
        run(genetic_distance_command, shell=True)
        print("--------------------------")


    def report_generation(self):
        """
        Report containing information about the identified sequences.
        """


    def multiqc(self, run_barcode):
        """
        Create MultiQC report, pulling in outputs from other tools
        """
        print("--------------------------\nMULTIQC\n--------------------------")
        if os.path.exists("{}/multiqc_report.html".format(self.multiqc_dir)):
            self.logger().warning("MULTIQC: Report already generated for {}".format(run_barcode))
        else:
            self.logger().info("MULTIQC: Creating MultiQC report for {}".format(run_barcode))
            multiqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} multiqc {} --outdir " \
                              "{}".format(self.app_dictionary["multiqc"], self.qc_dir, self.multiqc_dir)
            self.run_command(multiqc_command, "MULTIQC")
        print("--------------------------")


    def logger(self):
        """Write log messages to logfile.
        Arguments:
        message (str)
            Details about the logged event.
        tool (str)
            Tool name. Used to search within the insight ops website.
        """

        logging.basicConfig(format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG,
                            handlers=[
                                logging.FileHandler(self.logfile),
                                logging.StreamHandler(sys.stdout)]
                            )
        return logging


    def create_directory(self, directory):
        """
        Create output directories
        """
        if os.path.exists(directory):
            self.logger().warning("Directory not created - already exists: {}".format(directory))
        else:
            os.mkdir(directory)
            self.logger().info("Output directory created: {}".format(directory))
        return directory


    def get_centrifuge_index(self):
        """
        Download and unzip centrifuge index if not already downloaded
        """
        if not os.listdir(self.index_dir):
            centrifuge_download = "wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz " \
                                  "-p {} && tar xzvf {}/p_compressed+h+v.tar.gz -C {}".format(self.index_dir,
                                                                                              self.index_dir,
                                                                                              self.index_dir)
            self.run_command(centrifuge_download, "CENTRIFUGE INDEX")


    def run_command(self, command, tool):
        """"
        Run the command as a subprocess using Popen, write the stdout and stderr to file using logging, and write to
        log dependent on exitcode
        """
        # mlst and centrifuge-kreport need to separate stdout and stderr
        process = Popen(command, stdout=PIPE, stderr=STDOUT, shell=True)

        # per-line write to logfile
        for line in process.stdout:
            # decode and strip any leading/trailing whitespace
            line = line.decode("utf-8").strip()
            # remove control characters from string
            line = re.compile(r'[\n\r\t]').sub(' ', line)
            self.logger().info('got line from subprocess: {}'.format(line))

        exitcode = process.wait()

        success_message = "{}: Commmand completed successfully without any errors.".format(tool)
        failure_message = "{}: Job failed. See above subprocess output for details.".format(tool)
        if exitcode == 0:
            self.logger().info(success_message)
        else:
            self.logger().error(failure_message)

    #
    # def log_subprocess_output(self, out, err):
    #     for line in iter(out.readline, b''): # b'\n'-separated lines
    #
    #     for line in iter(err.readline, b''): # b'\n'-separated lines
    #         self.logger().error('got line from subprocess: %r', line)

def arg_parse():
    """
    Parses arguments supplied by the command line.
        :return: (Namespace object) parsed command line attributes
    Creates argument parser, defines command line arguments, then parses supplied command line arguments using the
    created argument parser.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--runfolder_dir', type=dir_path,
                        help="Directory containing runfolders to include in the analysis")
    parser.add_argument('-r', '--ref_dir', type=dir_path, help="Directory containing reference sequences")
    parser.add_argument('-i', '--centrifuge_index', type=dir_path, help="Directory containing centrifuge index")
    parser.add_argument('-o', '--out_dir', type=dir_path, help="Output directory")
    parser.add_argument('-g', '--ref_genome', type=file_path, help="Human reference genome fna file")

    return vars(parser.parse_args())

    # TEST
    # python3 pipeline.py
    # -d /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/test_run_folders
    # -r /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/reference_sequences
    # -i /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/centrifuge_index
    # -o /media/data2/share/outbreak_pipeline/rduffin/msc_project/test_output
    # -g /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna

    # LIVE
    # python3 pipeline.py
    # -d /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/run_folders
    # -r /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/reference_sequences
    # -i /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/centrifuge_index
    # -o /media/data2/share/outbreak_pipeline/rduffin/msc_project/output
    # -g /media/data2/share/outbreak_pipeline/rduffin/msc_project/data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna

def file_path(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_file:{path} is not a valid filepath")

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def main():
    args = arg_parse()
    runfolder_dir = args['runfolder_dir']
    reference_dir = args['ref_dir']
    index_dir = args['centrifuge_index']
    out_dir = args['out_dir']
    ref_gen = args['ref_genome']

    # Install containers
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])

    for run in os.listdir(runfolder_dir):
        # run the pipeline per run in runfolders directory
        analysis = Pipeline(run, runfolder_dir, reference_dir, index_dir, out_dir, ref_gen)
        analysis.set_off_analysis()

if __name__ == '__main__':
    main()

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
