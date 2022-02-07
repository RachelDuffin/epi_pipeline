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

class Pipeline:
    """Loop through and process runfolders in a given directory.
    A single class instance is required to process all runfolders.

    Methods defined here:
    loop_through_runs()

    """

    def __init__(self, run):
        """
        The constructor for Pipeline class
        """
        self.now = ""  # Stores time stamp for class instance, used in log file.
        self.base_dir = config.base_dir

        # specify data inputs paths
        self.reference_genome = config.reference_genome
        self.app_dictionary = app_dictionary
        self.run = run
        self.runfolder = "{}/{}".format(config.runfolders, self.run)
        self.sequencing_summary_file = glob.glob("{}/final_summary_*.txt".format(self.runfolder))
        self.run_id = self.get_identifier(file_list=self.sequencing_summary_file, string="sample_id=").rstrip("\n")
        self.fastq_path = "{}/*.fq".format(self.runfolder)
        self.reference_dir = "{}/reference_sequences/".format(self.base_dir)

        # specify outputs
        # self.out_dir = "{}/output/{}".format(os.getcwd(), self.run)
        self.out_dir = "{}/test_output/{}".format(os.getcwd(), self.run)
        self.centrifuge_index_summary = "{}/centrifuge_index_summary.txt".format(self.out_dir)
        self.out_subdirs = ["fastqc", "pycoqc", "human_read_removal", "de_novo_assembly", "mlst", "snp-sites",
                           "snp-dists", "centrifuge", "split_assemblies"]

        self.pycoqc_dir = "{}/pycoqc".format(self.out_dir)
        self.fastqc_dir = "{}/fastqc".format(self.out_dir)
        self.hr_rem_dir = "{}/human_read_removal".format(self.out_dir)
        self.class_dir = "{}/centrifuge".format(self.out_dir)

        self.index_dir = "{}/centrifuge_index".format(self.base_dir)
        self.index = "{}//p_compressed+h+v".format(self.index_dir)
        self.assembly_dir = "{}/de_novo_assembly".format(self.out_dir)
        self.assembly_split_dir = "{}/split_assemblies".format(self.out_dir)
        self.mlst_dir = "{}/mlst".format(self.out_dir)

        self.logfile = "{}/{}_logfile.txt".format(self.out_dir, self.run_id)


    def set_off_analysis(self):
        """
        Per-run analysis
        """
        if os.path.isdir(self.runfolder):
            # ensure log file is empty before starting processing
            open(self.logfile, "w")

            # create output directory for run
            self.create_directory(self.out_dir)

            # create subdirectories for outputs from each tool
            for subdir in self.out_subdirs:
                out_subdir = "{}/{}".format(self.out_dir, subdir)
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

                 if os.path.exists(assembly):
                     assembly_dict = self.assembly_split(run_barcode, assembly, assembly_outdir)

                     for contig in assembly_dict:
                         assembly_file = assembly_dict[contig]
                         run_barcode_contig = "{}_{}".format(run_barcode, contig)
                         self.mlst(run_barcode_contig, assembly_file)

                    # genetic_distance(out_dir=out_dir, run_id=run_id)
                    # variant_calling(out_dir=out_dir, run_id=run_id)
                    # multiqc(out_dir=out_dir, run_id=run_id)
                    # report_generation()

    def get_identifier(self, file_list, string):
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
        if not os.listdir(self.pycoqc_dir):
            # specify output directories
            self.split_barcodes()
            for file in glob.glob("{}/pycoqc/sequencing_summary_*".format(self.out_dir)):
                # get barcode name from barcode_split output file names
                barcode = file.split(".", 1)[0].split("summary_", 1)[1]
                self.logger().info("PYCOQC - Creating PycoQC json report for {}. ".format(self.run))
                self.logger().info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
                pycoqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                                 "{}/pycoqc/{}_pycoQC_output.json".format(self.app_dictionary["pycoqc"], file,
                                                                          self.out_dir, barcode)
                self.run_command(pycoqc_command, "PYCOQC")
        else:
            self.logger().info("PYCOQC - Directory not empty - pycoQC already run for {}".format(self.run))


    def split_barcodes(self):
        """
        Split barcodes
        """
        for file in glob.glob("{}/sequencing_summary_*".format(self.runfolder)):
            self.logger().info("PYCOQC BARCODE SPLIT: Splitting summary sequencing file {} "
                               "according to barcodes".format(file))
            self.logger().info("PYCOQC - version {}".format(self.app_dictionary["pycoqc"]))
            split_barcode_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
                                    "--output_unclassified --min_barcode_percent 0.0 --summary_file {} --output_dir " \
                                    "{}/pycoqc".format(self.app_dictionary["pycoqc"], file, self.out_dir)
            self.run_command(split_barcode_command, "PYCOQC BARCODE SPLIT")


    def fastqc(self, fastq, run_barcode):
        """
        FastQC analysis per barcode
        """
        fastqc_output = run_barcode + "_fastqc.html"
        # if output file already present, do not re-analyse
        if os.path.isfile("{}/{}".format(self.fastqc_dir, fastqc_output)):
            self.logger().warning("FASTQC: Output for {} already exists".format(run_barcode))
        else:
            self.logger().info("FASTQC: Creating FastQC file for {}".format(run_barcode))
            self.logger().info("FASTQC - version {}".format(self.app_dictionary["fastqc"]))
            fastqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                             "fastqc {} -o {}".format(self.app_dictionary["fastqc"], fastq, self.fastqc_dir)
            self.run_command(fastqc_command, "FASTQC")


    def human_read_removal(self, fastq, run_barcode):
        """
        Removes human reads from the samples by alignment to the human reference genome.
        """
        hr_rem_outpath = "{}/{}".format(self.hr_rem_dir, run_barcode)
        unmapped_fastq = "{}_unmapped.fastq".format(hr_rem_outpath)

        if os.path.isfile("{}_aligned.sam".format(hr_rem_outpath)):
            self.logger().warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            # align reads to human reference genome using ont-specific parameters
            self.logger().info("HR REMOVAL: Aligning reads to human reference genome for {}".format(fastq))
            self.logger().info("MINIMAP2 - version {}".format(self.app_dictionary["minimap2"]))
            minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} minimap2 -ax map-ont {} {} -o " \
                               "{}_aligned.sam".format(self.app_dictionary["minimap2"], self.reference_genome,
                                                       fastq, hr_rem_outpath)
            self.run_command(minimap2_command, "MINIMAP2")

        if os.path.isfile(unmapped_fastq):
            self.logger().warning("HR REMOVAL: Already complete for {}".format(run_barcode))
        else:
            # convert unassigned (non-human) reads from sam file to fastq file
            self.logger().info("SAMTOOLS FASTQ: Convert non-human from sam "
                               "to fastq for {}_aligned.sam".format(hr_rem_outpath))
            self.logger().info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))
            samtools_fastq_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} samtools fastq -f 4 " \
                                     "{}_aligned.sam -0 {}_unmapped.fastq".format(self.app_dictionary["samtools"],
                                                                                  hr_rem_outpath,
                                                                                  hr_rem_outpath)
            self.run_command(samtools_fastq_command, "SAMTOOLS FASTQ")

        if os.path.isfile("{}_samtools_stats.txt".format(hr_rem_outpath)):
            self.logger().warning("SAMTOOLS STATS: Samtools stats already conducted for {}".format(run_barcode))
                # remove intermediary file
                # os.remove("{}_unmapped.bam".format(out_path))
                # calculate % aligned reads to human reference genome
        else:
            self.logger().info("SAMTOOLS STATS: Run samtools stats for {}".format(run_barcode))
            self.logger().info("SAMTOOLS - version {}".format(self.app_dictionary["samtools"]))
            samtools_stats_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} samtools stats " \
                                     "{}_aligned.sam | grep ^SN | cut -f 2- > " \
                                     "{}_samtools_stats.txt".format(self.app_dictionary["samtools"],
                                                                    hr_rem_outpath,
                                                                    hr_rem_outpath)
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
            self.logger().info("CENTRIFUGE: Creating centrifuge index summary file for {}").format(run_barcode)
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

        if not os.path.isfile(centrifuge_kreport):
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

        if os.path.isfile(assembly_outfile_renamed):
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


    def assembly_split(self, run_barcode, assembly, assembly_outdir):
        """
        Split assembly output into a file per assembly using faidx.
        ADD FUNCTIONALITY TO ONLY KEEP CONTIGS THAT HAVE COVERAGE OVER
        """
        if os.listdir(self.assembly_split_dir):
            self.logger().warning("FAIDX: Assembly output already processed for {}".format(run_barcode))

        else:
            self.logger().info("FAIDX: Processing assembly output for {}".format(run_barcode))
            self.logger().info("PYFAIDX - version {}".format(self.app_dictionary["pyfaidx"]))
            faidx_command = "docker run -v `pwd`:`pwd` -v {}:{} {} /bin/sh -c 'cd {} ;" \
                            " faidx -x {}'".format(assembly_outdir, assembly_outdir, self.app_dictionary["pyfaidx"],
                                                  self.assembly_split_dir, assembly)
            self.run_command(faidx_command, "FAIDX")

        # ADD IN CODE TO PARSE THE assembly_info.txt FILE AND ONLY GET LINES THAT MATCH THE CONTIGS FROM THE TEXT FILE
        # WITH COVERAGE OVER 60X

        # create dictionary to store contig names and file paths
        assembly_dict = {}

        for file in os.listdir(self.assembly_split_dir):
            contig = file.split(".fasta")[0].rstrip("\n")
            filepath = "{}/{}".format(self.assembly_split_dir, file)
            assembly_dict[contig] = filepath

        return assembly_dict


    def mlst(self, run_barcode_contig, assembly_file):
        """
        Multi-locus sequence typing using MLSTcheck.
        Databases are bundled with the container so this doesn't require a separate download.
        """
        mlst_out_dir = "{}/{}".format(self.mlst_dir, run_barcode_contig)
        mlst_csv_output = "{}/mlst/{}.csv".format(self.out_dir, run_barcode_contig)

        directory = self.create_directory(mlst_out_dir)
        print("Directory: " + directory)
        print("Assembly file: " + assembly_file)

        if os.path.exists(mlst_csv_output):
            self.logger().warning("MLSTCHECK: mlst already complete for {}".format(run_barcode_contig))

        else:
            self.logger().info("MLSTCHECK: Conducting MLST for {}".format(run_barcode_contig))
            self.logger().info("MLSTCHECK - version {}".format(self.app_dictionary["mlstcheck"]))

            # no -s option provided to get_sequence_type meaning every databases within the database snapshot in the docker
            # will be searched
            mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -it {} get_sequence_type {} -o " \
                           "{}".format(self.app_dictionary["mlstcheck"], assembly_file, directory)
            print(mlst_command)
            # self.run_command(mlst_command, "MLSTCHECK")
        # docker run --rm -v `pwd`:`pwd` -w `pwd` -it quay.io/biocontainers/perl-bio-mlst-check@sha256:8
        # 28f4a536603a39559118a279e9c66d71e83a155cf0ac824efdb9339ba59e201 get_sequence_type
        # /media/data2/share/outbreak_pipeline/rduffin/msc_project/output/200408_CMG_RUN3_Pool3/split_assemblies/200408_CMG_RUN3_Pool3_barcode01_contig_10_processed.fasta -o
        # /media/data2/share/outbreak_pipeline/rduffin/msc_project/output/200408_CMG_RUN3_Pool3/mlst/200408_CMG_RUN3_Pool3_barcode01_contig_10


    def variant_calling(self):
        """
        Identify SNPs between samples
        """
        print("--------------------------\nVARIANT CALLING\n--------------------------")
        print(self.out_dir)
        print(self.run_id)
        # SNP-sites to identify SNPs between samples
        variant_calling_command = "docker run --rm -v `pwd`:`pwd` -w `{} snp-sites -m -o {} " \
                                  "{}".format(self.app_dictionary["snp-sites"], OUTPUT_FILENAME, INPUT_FILENAME)
        run(variant_calling_command, shell=True)
        print("--------------------------")
        pass


    def genetic_distance(self):
        """
        Calculate SNP distances
        """
        print("--------------------------\nGENETIC DISTANCE CALCULATION\n--------------------------")
        # SNP-dists to calculate SNP distances
        genetic_distance_command = "docker run --rm -v `pwd`:`pwd` -w `{} snp-dists {} > " \
                                   "{}".format(self.app_dictionary["snp-dists"], INPUT_FILENAME, OUTPUT_FILENAME.tsv)
        run(genetic_distance_command, shell=True)
        print("--------------------------")
        pass


    def report_generation(self):
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
                              "{}/multiqc/{}".format(self.app_dictionary["multiqc"], out_dir, out_dir, run_id)
            run(multiqc_command, shell=True)
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
        #
        #
        # logging.basicCo
        # nfig(filename=self.logfile,
        #                     filemode='a',
        #                     format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
        #                     datefmt='%H:%M:%S',
        #                     level=logging.DEBUG,
        #                     )
        return logging


    def create_directory(self, directory):
        """
        Create output directories
        """
        if os.path.exists(directory):
            self.logger().warning("Directory not created - already exists: {}".format(directory))
        else:
            print(directory)
            os.mkdir(directory)
            self.logger().info("Output directory created: {}".format(directory))
        return directory


    def get_centrifuge_index(self):
        """
        Download and unzip centrifuge index if not already downloaded
        """
        if not os.listdir(self.index_dir):
            centrifuge_download = "wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz -p " \
                                  "{} && tar xzvf {}/p_compressed+h+v.tar.gz -C {}".format(self.index_dir, self.index_dir,
                                                                                           self.index_dir)
            self.run_command(centrifuge_download, "CENTRIFUGE INDEX")


    def run_command(self, command, tool):
        """"
        Run the command as a subprocess using Popen, write the stdout and stderr to file using logging, and write to
        log dependent on exitcode
        """
        process = Popen(command, stdout=PIPE, stderr=STDOUT, shell=True)

        # per-line write to logfile
        for line in process.stdout:
            # decode and strip any leading/trailing whitespace
            line = line.decode("utf-8").strip()
            # remove control characters from string
            line = re.compile(r'[\n\r\t]').sub(' ', line)
            print([line])
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

def main():

    # Install containers
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])
        
    for run in os.listdir(config.runfolders):
        # run the pipeline per run in runfolders directory
        analysis = Pipeline(run)
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
