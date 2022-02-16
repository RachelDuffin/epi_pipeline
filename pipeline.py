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
        self.out_subdirs = ["QC", "human_read_removal", "de_novo_assembly", "split_assemblies", "snp-sites",
                           "snp-dists", "centrifuge", "split_barcodes", "passing_contigs"]

        self.hr_rem_dir = "{}/human_read_removal".format(self.runfolder_out)
        self.class_dir = "{}/centrifuge".format(self.runfolder_out)

        self.index = "{}//p_compressed+h+v".format(self.index_dir)
        self.assembly_dir = "{}/de_novo_assembly".format(self.runfolder_out)
        self.assembly_split_dir = "{}/split_assemblies".format(self.runfolder_out)
        self.qc_dir = "{}/QC".format(self.runfolder_out)
        self.split_barcodes_dir = "{}/split_barcodes".format(self.runfolder_out)
        self.passing_contigs_dir = "{}/passing_contigs".format(self.runfolder_out)

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

            logger(self.logfile).info("ANALYSING RUN {}".format(self.run))
            # split barcodes for whole run
            self.split_barcodes()

            for fastq in glob.glob(self.fastq_path):

                barcode = self.get_identifier(file_list=[fastq], string="barcode")
                run_barcode = "{}_{}".format(self.run_id, barcode)

                # run pre-alignment QC tools
                self.pycoqc(run_barcode, barcode)
                self.fastqc(fastq, run_barcode, barcode)

                # conduct human read removal
                unmapped_fastq = self.human_read_removal(fastq, run_barcode)
                self.classification(run_barcode, unmapped_fastq)
                assembly, assembly_outdir = self.de_novo_assembly(run_barcode, unmapped_fastq)

                # get dictionary of contigs that pass coverage requirements
                passing_contigs = self.get_passing_contigs(run_barcode, assembly_outdir)

                # get list of all filepaths containing passing contigs
                assembly_list = self.assembly_split(run_barcode, assembly, assembly_outdir, passing_contigs)
                # concatenate into a single passing contigs file
                passing_contigs_fq = self.concatenate_sequences(assembly_list, run_barcode)

                fq_headers_renamed = self.rename_headers(run_barcode, passing_contigs_fq)

                    # genetic_distance(out_dir=out_dir, run_id=run_id)
                    # variant_calling(out_dir=out_dir, run_id=run_id)
                    # report_generation()
        return


    def get_passing_contigs(self, run_barcode, assembly_outdir):
        """
        Parse assembly info file to get a list of contigs that pass minimum coverage requirements (60X)
        """
        contig_list = []
        assembly_info = "{}/assembly_info.txt".format(assembly_outdir)

        if os.path.isfile(assembly_info):
            logger(self.logfile).info("PARSE ASSEMBLY INFO: Parsing assembly info file for {}".format(run_barcode))
            with open(assembly_info, "r") as file:
                for line in file.readlines():
                    if not line.startswith("#seq_name"):
                        line = line.split("\t")
                        coverage = int(line[2])
                        if coverage >= 60:
                            contig_list.append("{}.".format(line[0]))
        else:
            logger(self.logfile).warning("PARSE ASSEMBLY INFO: Assembly info "
                                         "file does not exist for {}".format(run_barcode))
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


    def pycoqc(self, run_barcode, barcode):
        """
        Then creates a pycoQC json report per barcode.
        """
        # If pycoQC not yet run, run pycoQC
        sequencing_summary_file = "{}/sequencing_summary_{}.txt".format(self.split_barcodes_dir, barcode)
        pycoqc_outfile = "{}/{}_pycoQC_output.json".format(self.qc_dir, run_barcode)

        if os.path.isfile(sequencing_summary_file):
            if not os.path.isfile(pycoqc_outfile):
                logger(self.logfile).info("PYCOQC - Creating PycoQC json report for {}. ".format(run_barcode))
                logger(self.logfile).info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
                pycoqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                                 "{}".format(self.app_dictionary["pycoqc"],
                                             sequencing_summary_file, pycoqc_outfile)
                run_command(self.logfile, pycoqc_command, "PYCOQC")
            else:
                logger(self.logfile).info("PYCOQC - Directory not empty - pycoQC "
                                          "already run for {}".format(run_barcode))
        else:
            logger(self.logfile).info("PYCOQC - Summary file for {} could not be found. ".format(run_barcode))


    def split_barcodes(self):
        """
        Split barcodes. Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes()
        function.
        """
        if not glob.glob("{}/sequencing_summary_*".format(self.split_barcodes_dir)):

            for file in glob.glob("{}/sequencing_summary_*".format(self.runfolder)):
                logger(self.logfile).info("PYCOQC BARCODE SPLIT: Splitting summary sequencing file {} "
                                   "according to barcodes".format(file))
                logger(self.logfile).info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
                split_barcode_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
                                        "--output_unclassified --min_barcode_percent 0.0 --summary_file {} " \
                                        "--output_dir {}".format(self.app_dictionary["pycoqc"], file,
                                                                 self.split_barcodes_dir)
                run_command(self.logfile, split_barcode_command, "PYCOQC BARCODE SPLIT")

        else:
            logger(self.logfile).warning("PYCOQC BARCODE SPLIT: Summary "
                                  "sequencing file already split for {} ".format(self.run_id))


    def fastqc(self, fastq, run_barcode, barcode):
        """
        FastQC analysis per barcode
        """
        # if output file already present, do not re-analyse
        if glob.glob("{}/*{}*.html".format(self.qc_dir, barcode)):
            logger(self.logfile).warning("FASTQC: Output for {} already exists".format(run_barcode))
        else:
            logger(self.logfile).info("FASTQC: Creating FastQC file for {}".format(run_barcode))
            logger(self.logfile).info("FASTQC - version {}".format(self.app_dictionary["fastqc"]))
            fastqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                             "fastqc {} -o {}".format(self.app_dictionary["fastqc"], fastq, self.qc_dir)
            run_command(self.logfile, fastqc_command, "FASTQC")


    def human_read_removal(self, fastq, run_barcode):
        """
        Removes human reads from the samples by alignment to the human reference genome, with output in sam format.
        Converts the sam output to
        Runs samtools stats on the alignment output to show how many reads mapped to the human genome.
        """
        hr_rem_outpath = "{}/{}".format(self.hr_rem_dir, run_barcode)
        unmapped_fastq = "{}_unmapped.fastq".format(hr_rem_outpath)
        sam_output = "{}_aligned.sam".format(hr_rem_outpath)
        samtools_stats_file = "{}/{}_samtools_stats.txt".format(self.qc_dir, run_barcode)

        # align reads to human reference genome using ont-specific parameters
        # -a generates the outputs in the SAM format which allows extraction of unmapped reads
        # -x -map-ont aligns long noisy reads of ~10% error rate to reference genome

        if os.path.isfile(sam_output):
            logger(self.logfile).warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            logger(self.logfile).info("HR REMOVAL: Aligning reads to human reference genome for {}".format(fastq))
            logger(self.logfile).info("MINIMAP2 - version {}".format(self.app_dictionary["minimap2"]))

            minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} minimap2 -ax " \
                               "map-ont {} {} -o {}".format(self.app_dictionary["minimap2"],
                                                            self.reference_genome, fastq, sam_output)
            run_command(self.logfile, minimap2_command, "MINIMAP2")

        # Perform samtools stats on the alignment output to calculate alignment stats
        # (eg. no. reads mapped to human genome)
        if os.path.isfile(samtools_stats_file):
            logger(self.logfile).warning("SAMTOOLS STATS: Samtools stats already conducted for {}".format(run_barcode))
        else:
            logger(self.logfile).info("SAMTOOLS STATS: Run samtools stats for {}".format(run_barcode))
            logger(self.logfile).info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))

            samtools_stats_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} samtools stats " \
                                     "{} > {}".format(self.app_dictionary["samtools"], sam_output, samtools_stats_file)
            run_command(self.logfile, samtools_stats_command, "SAMTOOLS STATS")

        # Convert sam file to fastq file using samtools, outputting only unmapped reads (non-human)
        # -f 4 = only outputs unmapped segments (with bitwise flag integer = 4)
        if os.path.isfile(unmapped_fastq):
            logger(self.logfile).warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            logger(self.logfile).info("SAMTOOLS FASTQ: Convert non-human from sam "
                                      "to fastq for {}_unmapped.sam".format(hr_rem_outpath))
            logger(self.logfile).info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))
            samtools_fastq_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} " \
                                     "samtools fastq -f 4 {} -0 {}".format(self.app_dictionary["samtools"],
                                                                            sam_output, unmapped_fastq)
            run_command(self.logfile, samtools_fastq_command, "SAMTOOLS FASTQ")

        return unmapped_fastq


    def concatenate_sequences(self, filepath_list, run_barcode):
        """
        Concatenate input sequences into one fna file
        """
        passing_contigs_fq = "{}/{}_passing_assembly.fasta".format(self.passing_contigs_dir, run_barcode)

        if filepath_list:
            if os.path.exists(passing_contigs_fq):
                logger(self.logfile).warning("CONCATENATE SEQUENCES: Passing sequencings already "
                                                  "concatenated for {}".format(run_barcode))
            else:
                logger(self.logfile).info("CONCATENATE SEQUENCES: Concatenating passing "
                                               "sequences for {}".format(run_barcode))
                assembly_string = ""
                for filepath in filepath_list:
                    assembly_string = "{} {}".format(assembly_string, filepath)

                concat_command = "cat {} > {}".format(assembly_string, passing_contigs_fq)
                run_command(self.logfile, concat_command, "CONCATENATE SEQUENCES")
        else:
            logger(self.logfile).warning("CONCATENATE SEQUENCES: No passing sequences found for {}".format(run_barcode))

        return passing_contigs_fq


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
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge index summary file already created")
        else:
            logger(self.logfile).info("CENTRIFUGE: Creating centrifuge index summary file for {}".format(run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            # run centrifuge-inspect to output a summary of the index used for classification
            centrifuge_inspect = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} centrifuge-inspect -s " \
                                 "{} > {}".format(self.app_dictionary["centrifuge"], self.index,
                                                  self.centrifuge_index_summary)
            run_command(self.logfile, centrifuge_inspect, "CENTRIFUGE")

        # Run centrifuge classification
        if os.path.isfile(centrifuge_results) and os.path.isfile(centrifuge_summary):
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge already run for {}".format(run_barcode))
        else:
            # classify reads
            logger(self.logfile).info("CENTRIFUGE: Running centrifuge classification for {}".format(run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            # -S is file to write classification results to, --report-file is file to write classification summary to,
            # -x is the index, -f is the input sequence, --env sets centriuge indexes environment variable so it knows
            # where to look for the index
            centrifuge_command = "docker run --env CENTRIFUGE_INDEXES={} --rm -v `pwd`:`pwd` -w `pwd` {} " \
                                 "centrifuge -x p_compressed+h+v -q {} -S {} " \
                                 "--report-file {}".format(self.index_dir, self.app_dictionary['centrifuge'],
                                                           unmapped_fastq, centrifuge_results, centrifuge_summary)
            run_command(self.logfile, centrifuge_command, "CENTRIFUGE")

        # Convert centrifuge report to a kreport (kraken-style report). This format is picked up by MultiQC, and can
        # be loaded into Pavian
        if os.path.isfile(centrifuge_kreport):
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge-kreport already generated for {}".format(run_barcode))
        else:
            logger(self.logfile).info("CENTRIFUGE: Generating centrifuge-kreport for {}".format(run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(self.app_dictionary["centrifuge"]))
            kreport_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} centrifuge-kreport -x {}" \
                              " {} > {}".format(self.app_dictionary["centrifuge"], self.index,
                                                centrifuge_results, centrifuge_kreport)
            run_command(self.logfile, kreport_command, "CENTRIFUGE")


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
            logger(self.logfile).warning("FLYE: De novo assembly already complete for {}".format(run_barcode))
        else:
            logger(self.logfile).info("FLYE: Performing FLYE de novo assembly for {}".format(unmapped_fastq))
            logger(self.logfile).info("FLYE - version {}".format(self.app_dictionary["flye"]))
            self.create_directory(assembly_outdir)
            flye_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} flye --nano-raw {} --out-dir {} " \
                           "--meta".format(self.app_dictionary["flye"], unmapped_fastq, assembly_outdir)

            # rename output to add in run ID and barcode
            run_command(self.logfile, flye_command, "FLYE")

            if os.path.exists(assembly_outfile):
                try:
                    os.rename(assembly_outfile, assembly_outfile_renamed)
                except:
                    logger(self.logfile).error("FLYE: Assembly outfile {} could "
                                               "not be renamed".format(assembly_outfile))
            else:
                logger(self.logfile).warning("FLYE: No assembly produced.")

        return assembly_outfile_renamed, assembly_outdir


    def assembly_split(self, run_barcode, assembly, assembly_outdir, passing_contigs):
        """
        Split assembly output into a file per assembly using faidx.
        ADD FUNCTIONALITY TO ONLY KEEP CONTIGS THAT HAVE COVERAGE OVER
        """
        if os.listdir(self.assembly_split_dir):
            logger(self.logfile).warning("FAIDX: Assembly output already processed for {}".format(run_barcode))

        else:
            if os.path.exists(assembly):
                logger(self.logfile).info("FAIDX: Processing assembly output for {}".format(run_barcode))
                logger(self.logfile).info("PYFAIDX - version {}".format(self.app_dictionary["pyfaidx"]))
                faidx_command = "docker run -v `pwd`:`pwd` -v {}:{} {} /bin/sh -c 'cd {} ; " \
                                "faidx -x {}'".format(assembly_outdir, assembly_outdir, self.app_dictionary["pyfaidx"],
                                                      self.assembly_split_dir, assembly)

                run_command(self.logfile, faidx_command, "FAIDX")
            else:
                logger(self.logfile).info("FAIDX: Assembly file not present for {}".format(run_barcode))

        # create dictionary to store only filepaths of contigs that surpass the coverage threshold of 60X
        assembly_list = []

        for file in os.listdir(self.assembly_split_dir):
            filepath = "{}/{}".format(self.assembly_split_dir, file)
            if any(contig in file for contig in passing_contigs):
                assembly_list.append(filepath)

        return assembly_list


    def rename_headers(self, run_barcode, passing_contigs_fq):
        """
        Add run and barcode name into the read headers in fasta files
        """
        newfasta = "{}_headers_renamed.fasta".format(passing_contigs_fq.split(".fasta")[0])
        if os.path.isfile(passing_contigs_fq):
            if os.path.isfile(newfasta):
                logger(self.logfile).warning("RENAME HEADERS: Headers already renamed for {}".format(run_barcode))

            else:
                logger(self.logfile).warning("RENAME HEADERS: Renaming headers for {}".format(run_barcode))

                fasta = open(passing_contigs_fq, "r")
                newfasta = open(newfasta, "w")

                for line in fasta.readlines():
                    if line.startswith(">"):
                        elements = line.split(">")
                        new_header = ">{}{}_{}".format(elements[0], run_barcode, elements[1])
                        newfasta.write(new_header)
                    else:
                        newfasta.write(line)
                fasta.close()
                newfasta.close()
        else:
            logger(self.logfile).warning("RENAME HEADERS: No passing sequences found for {}".format(run_barcode))


        return newfasta


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




    def report_generation(self):
        """
        Report containing information about the identified sequences.
        """


    def create_directory(self, directory):
        """
        Create output directories
        """
        if os.path.exists(directory):
            logger(self.logfile).warning("Directory not created - already exists: {}".format(directory))
        else:
            os.mkdir(directory)
            logger(self.logfile).info("Output directory created: {}".format(directory))
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
            run_command(self.logfile, centrifuge_download, "CENTRIFUGE INDEX")

    #
    # def log_subprocess_output(self, out, err):
    #     for line in iter(out.readline, b''): # b'\n'-separated lines
    #
    #     for line in iter(err.readline, b''): # b'\n'-separated lines
    #         logger(self.logfile).error('got line from subprocess: %r', line)

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


def file_path(path):
    """
    Checks the command line argument provided is an existing filepath.
    """
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_file:{path} is not a valid filepath")


def dir_path(path):
    """
    Checks the command line argument provided is an existing directory.
    """
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def logger(logfile):
    """Write log messages to logfile.
    Arguments:
    message (str)
        Details about the logged event.
    tool (str)
        Tool name. Used to search within the insight ops website.
    """
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                            logging.FileHandler(logfile),
                            logging.StreamHandler(sys.stdout)]
                        )
    return logging


def mlst(out_dir, assembly_files, logfile):
    """
    Multi-locus sequence typing using mlst (tseemann/mlst).
    Databases are bundled with the container so this doesn't require a separate download.
    """
    mlst_output = "{}/mlst.tsv".format(out_dir)
    mlst_scheme_list = "{}/mlst_schemes.txt".format(out_dir)

    if os.path.exists(mlst_output):
        logger(logfile).warning("MLST: mlst already complete (across all runs)")
    else:
        # Record pubmlst schemes
        logger(logfile).info("MLST - recording supported PubMLST schemes in file")

        schemes_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} mlst " \
                          "--longlist > {}".format(app_dictionary["mlst"], mlst_scheme_list)

        run_command(logfile, schemes_command, "MLST")

        # Run mlst across all runs
        logger(logfile).info("MLST: Conducting MLST across all runs")
        logger(logfile).info("MLST - version {}".format(app_dictionary["mlst"]))

        mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {}" \
                       " mlst {} > {}".format(app_dictionary["mlst"], assembly_files, mlst_output)

        run_command(logfile, mlst_command, "MLST")


def multiqc(out_dir, logfile):
    """
    Create MultiQC report, pulling in outputs from other tools
    """
    if os.path.exists("{}/multiqc_report.html".format(out_dir)):
        logger(logfile).warning("MULTIQC: Report already generated")
    else:
        logger(logfile).info("MULTIQC: Creating MultiQC report")
        multiqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} multiqc {} --outdir " \
                          "{}".format(app_dictionary["multiqc"], out_dir, out_dir)
        run_command(logfile, multiqc_command, "MULTIQC")


def run_command(logfile, command, tool):
    """"
    Run the command as a subprocess using Popen, write the stdout and stderr to file using logging, and write to
    log dependent on exitcode
    """
    process = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    stdout, stderr = process.communicate()
    exitcode = process.wait()

    success_message = "{}: Command completed successfully without any errors.".format(tool)
    failure_message = "{}: Job failed. See above subprocess output for details.".format(tool)

    if exitcode == 0:
        logger(logfile).info(success_message)
    else:
        logger(logfile).error(failure_message)

    for stream in [stdout, stderr]:
        # decode and split by line
        stream = stream.decode("utf-8")
        lines = stream.split("\n")
        for line in lines:
        #     # remove control characters from string
        #     line = re.compile(r'[\n\r\t]').sub(' ', line)
            logger(logfile).info('got line from subprocess: {}'.format(line))


def multi_fasta_alignment(out_dir, assembly_files, analysis_logfile):
    """
    Create multi-fasta alignment across all sequences using ??? minimap2
    """
    concatenated_assemblies = "{}/concatenated_assemblies.fasta".format(out_dir)

    alignment_multifasta = "{}/multi_fasta_alignment.fastq".format(out_dir)

    # concatenate contents of assembly files
    if os.path.isfile(concatenated_assemblies):
        logger(analysis_logfile).warning("CONCATENATE SEQUENCES: Already concatenated sequences across the analysis")

    else:
        logger(analysis_logfile).info("CONCATENATE SEQUENCES: Concatenating sequences "
                                      "into a single fasta file across the analysis")
        assembly_string = ""

        for filepath in glob.glob(assembly_files):
            assembly_string = "{} {}".format(assembly_string, filepath)

        concat_command = "cat {} > {}".format(assembly_string, concatenated_assemblies)
        run_command(analysis_logfile, concat_command, "CONCATENATE SEQUENCES")

    # alignment_sam = "{}/multi_fasta_alignment.sam".format(out_dir)


    # creating multi fasta alignment from all filtered alignment files
    if os.path.isfile(alignment_sam):
        logger(analysis_logfile).warning("MULTI FASTA ALIGNMENT: Already complete (across all runs)")
    else:
        logger(analysis_logfile).info("MULTI FASTA ALIGNMENT: Performing all-vs-all alignment across all runs")
        logger(analysis_logfile).info("MINIMAP2 - version {}".format(app_dictionary["minimap2"]))

        minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} minimap2 -ax " \
                           "map-ont {} -o {}".format(app_dictionary["minimap2"], assembly_files, alignment_sam)

        run_command(analysis_logfile, minimap2_command, "MINIMAP2")

    # converting sam multi alignment file to fastq
    if os.path.isfile(alignment_multifasta):
        logger(analysis_logfile).warning("MULTI FASTA ALIGNMENT: Sam already converted to fastq")
    else:
        logger(analysis_logfile).info("MULTI FASTA ALIGNMENT: Converting sam all vs. all alignment to fastq")

        samtools_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} samtools " \
                           "bam2fq {} > {}".format(app_dictionary["samtools"], alignment_sam, alignment_multifasta)

        run_command(analysis_logfile, samtools_command, "SAMTOOLS")


def genetic_distance(out_dir, analysis_logfile):
    """
    Calculate SNP distances from a multi-fasta alignment file.
    """
    multi_alignment = "{}/multi_fasta_alignment.sam".format(out_dir)
    snp_matrix = "{}/snp_matrix.tsv".format(out_dir)

    if os.path.isfile(snp_matrix):
        logger(analysis_logfile).warning("SNP DISTANCE CALCULATION: Already complete (across all runs)")
    else:
        logger(analysis_logfile).info("SNP DISTANCE CALCULATION: Performing snp distance "
                                      "calculation on multi fasta alignment")
        logger(analysis_logfile).info("SNP-DISTS - version {}".format(app_dictionary["snp-dists"]))

    # GENETIC DISTANCE CALCULATION
    # SNP-dists to calculate SNP distances
    genetic_distance_cmd = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} snp-dists {} > " \
                               "{}".format(app_dictionary["snp-dists"], multi_alignment, snp_matrix)

    run_command(analysis_logfile, genetic_distance_cmd, "SNP-DISTS")


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

    # Conduct MLST across all runfolders in the analysis
    analysis_logfile = "{}/analysis_logfile.txt".format(out_dir)
    assembly_files = "{}/*/passing_contigs/*passing_assembly_headers_renamed.fasta".format(out_dir)

    mlst(out_dir, assembly_files, analysis_logfile)

    multi_fasta_alignment(out_dir, assembly_files, analysis_logfile)


    genetic_distance(out_dir, analysis_logfile)

    # Collect QC data for multiqc for a single MultiQC report for all samples in the analysis
    qc_dirs = "{}/*/QC/".format(out_dir)
    multiqc(out_dir, analysis_logfile)


if __name__ == '__main__':
    main()