"""
Epidemiology pipeline

Pipeline is structured across two main classes:

    One at the analysis level (across all runfolders) which run tools
        that produce a single output per analysis e.g. sequence typing and multiqc report generation.
    One at the runfolder level which loops through samples at the runfolder level. e.g. human read removal,
        classification, de novo assembly.

Instructions on how to use this script can be found in the Readme file.
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


class AllRuns:
    """
    Conduct processing across all runfolders.

    Methods are defined here:
        set_off_analysis()          Calls analysis steps to be conducted across all runfolders,
                                    then sets of RunfolderAnalysis (per-runfolder analysis)
        get_centrifuge_index()      Download and unzip centrifuge index if not already downloaded
        mlst()                      Multi-locus sequence typing using mlst (tseemann/mlst).
        multiqc()                   Create MultiQC report for the analysis
    """

    def __init__(self, runfolder_dir, reference_dir, index_dir, out_dir, ref_gen):
        """
        The constructor for AllRuns class
        """
        self.analysis_logfile = "{}/analysis_logfile.txt".format(out_dir)
        self.runfolder_dir = runfolder_dir
        self.reference_dir = reference_dir
        self.index_dir = index_dir
        self.out_dir = out_dir
        self.ref_gen = ref_gen
        self.refseq_dict = config.refseq_dict
        self.assembly_files_filtered = "{}/*/passing_contigs/*passing_assembly_headers_renamed.fasta".format(
            self.out_dir)
        self.assembly_files_unfiltered = "{}/*/de_novo_assembly/*/*_assembly.fasta".format(self.out_dir)
        self.filtered_mlst_output = "{}/coverage_filtered_mlst.tsv".format(self.out_dir)
        self.unfiltered_mlst_output = "{}/unfiltered_mlst.tsv".format(out_dir)
        self.mlst_scheme_list = "{}/mlst_schemes.txt".format(self.out_dir)
        self.multiqc_search_paths_string = " ".join(['{}/*/QC/'.format(self.out_dir),
                                                     '{}/*/centrifuge/'.format(self.out_dir)])
        self.multiqc_config = '{}/multiqc_config.yaml'.format(os.getcwd())


    def set_off_analysis(self):
        """
        Analysis across all supplied runfolders
        """
        # Install containers
        for key in app_dictionary:
            install_containers.install_tools(key, app_dictionary[key])

        # get any reference material that isn't already present
        self.get_centrifuge_index()
        get_reference_sequences.download_sequences(self.refseq_dict, self.reference_dir)

        for run in os.listdir(self.runfolder_dir):
            # run the pipeline per run in runfolders directory
            runfolder_analysis = RunfolderAnalysis(run, self.runfolder_dir, self.index_dir, self.out_dir, self.ref_gen)
            runfolder_analysis.set_off_analysis()

        # Conduct MLST across all runfolders in the analysis, on both coverage-filtered and unfiltered outputs
        self.mlst(self.assembly_files_filtered, self.filtered_mlst_output)
        # conduct mlst on non-coverage-filtered assemblies
        self.mlst(self.assembly_files_unfiltered, self.unfiltered_mlst_output)

        self.multiqc()


    def get_centrifuge_index(self):
        """
        Download and unzip centrifuge index if not already downloaded
        """
        if not os.listdir(self.index_dir):
            centrifuge_download = "wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz " \
                                  "-p {} && tar xzvf {}/p_compressed+h+v.tar.gz -C {}".format(self.index_dir,
                                                                                              self.index_dir,
                                                                                              self.index_dir)
            run_command(self.analysis_logfile, centrifuge_download, "CENTRIFUGE INDEX")


    def mlst(self, assembly_files, mlst_output):
        """
        Multi-locus sequence typing using mlst (tseemann/mlst).
        Databases are bundled with the container so this doesn't require a separate download.
        """
        if os.path.exists(mlst_output):
            logger(self.analysis_logfile).warning("MLST: mlst already complete (across all runs)")
        else:
            # Record pubmlst schemes
            logger(self.analysis_logfile).info("MLST - recording supported PubMLST schemes in file")

            schemes_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} mlst " \
                              "--longlist > {}".format(app_dictionary["mlst"], self.mlst_scheme_list)

            run_command(self.analysis_logfile, schemes_command, "MLST")

            # Run mlst across all runs
            logger(self.analysis_logfile).info("MLST: Conducting MLST across all runs")
            logger(self.analysis_logfile).info("MLST - version {}".format(app_dictionary["mlst"]))

            mlst_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {}" \
                           " mlst {} > {}".format(app_dictionary["mlst"], assembly_files, mlst_output)

            run_command(self.analysis_logfile, mlst_command, "MLST")

    def multiqc(self):
        """
        Create MultiQC report, pulling in outputs from other tools
        """
        if os.path.exists("{}/multiqc_report.html".format(self.out_dir)):
            logger(self.analysis_logfile).warning("MULTIQC: Report already generated")
        else:
            logger(self.analysis_logfile).info("MULTIQC: Creating MultiQC report")
            multiqc_command = "docker run --rm -v `pwd`:`pwd` -v {}:`pwd`/multiqc_config.yaml -w `pwd` {} multiqc {} " \
                              "--outdir {}".format(self.multiqc_config, app_dictionary["multiqc"],
                                                   self.multiqc_search_paths_string, self.out_dir)
            run_command(self.analysis_logfile, multiqc_command, "MULTIQC")


class RunfolderAnalysis:
    """
    Loop through and process runfolders in a given directory.
    A class instance is required per runfolder.

    Methods defined here:
        set_off_analysis()          Calls per-run analysis steps, then sets off FastqAnalysis (per-fastq)
        split_barcodes()            Splits summary sequence file according to barcodes
        get_identifier()            Get run id or barcode from sequencing summary file or fastq file
    """

    def __init__(self, run, runfolder_dir, index_dir, out_dir, ref_gen):
        """
        The constructor for RunfolderAnalysis class
        """
        # specify data inputs paths
        self.index_dir = index_dir
        self.reference_genome = ref_gen
        self.runfolder_dir = runfolder_dir
        self.run = run
        self.runfolder = "{}/{}".format(self.runfolder_dir, self.run)
        self.fastq_path = "{}/*.fq".format(self.runfolder)
        self.final_summary_file = glob.glob("{}/final_summary_*.txt".format(self.runfolder))
        self.run_id = self.get_identifier(self.final_summary_file, "sample_id=").rstrip("\n")

        # specify output paths
        self.out_dir = out_dir
        self.runfolder_out = "{}/{}".format(out_dir, self.run)
        self.out_subdirs = ["QC", "human_read_removal", "de_novo_assembly", "split_assemblies", "snp-sites",
                           "snp-dists", "centrifuge", "split_barcodes", "passing_contigs"]
        self.logfile = "{}/{}_logfile.txt".format(self.runfolder_out, self.run_id)
        self.split_barcodes_dir = "{}/split_barcodes".format(self.runfolder_out)
    
    
    def set_off_analysis(self):
        """
        Per-run analysis
        """
        if os.path.isdir(self.runfolder):
            # create run output dir if it doesn't already exist
            create_directory(self.runfolder_out, self.logfile)

            # create logfile if it doesn't already exist
            if not os.path.isfile(self.logfile):
                open(self.logfile, 'w')

            # create output directory for run
            create_directory(self.out_dir, self.logfile)

            # create subdirectories for outputs from each tool
            for subdir in self.out_subdirs:
                out_subdir = "{}/{}".format(self.runfolder_out, subdir)
                create_directory(out_subdir, self.logfile)

            logger(self.logfile).info("ANALYSING RUN {}".format(self.run))
            # split barcodes for whole run, so that pycoQC can be conducted only on the relevant samples
            self.split_barcodes()

            for fastq in glob.glob(self.fastq_path):

                barcode = self.get_identifier(file_list=[fastq], string="barcode")
                run_barcode = "{}_{}".format(self.run_id, barcode)

                fastq_analysis = FastqAnalysis(fastq, self.runfolder, self.run_id, barcode, run_barcode, self.index_dir,
                                               self.split_barcodes_dir, self.runfolder_out, self.logfile,
                                               self.reference_genome)
                fastq_analysis.set_off_analysis()


    def split_barcodes(self):
        """
        Split barcodes. Using guppy barcoding file. Checks if barcodes already split, if not calls split_barcodes()
        function.
        """
        if not glob.glob("{}/sequencing_summary_*".format(self.split_barcodes_dir)):

            for file in glob.glob("{}/sequencing_summary_*".format(self.runfolder)):
                logger(self.logfile).info("PYCOQC BARCODE SPLIT: Splitting summary sequencing file {} "
                                   "according to barcodes".format(file))
                logger(self.logfile).info("PYCOQC - version {}".format(app_dictionary["pycoqc"]))
                split_barcode_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} Barcode_split " \
                                        "--output_unclassified --min_barcode_percent 0.0 --summary_file {} " \
                                        "--output_dir {}".format(app_dictionary["pycoqc"], file,
                                                                 self.split_barcodes_dir)
                run_command(self.logfile, split_barcode_command, "PYCOQC BARCODE SPLIT")

        else:
            logger(self.logfile).warning("PYCOQC BARCODE SPLIT: Summary "
                                         "sequencing file already split for {} ".format(self.run_id))

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


class FastqAnalysis:
    """
    Loop through and process fastqs in a runfolder.
    A class instance is required per fastq.

    Methods defined here:
        pycoqc()                    Create pycoQC json report for barcode (sample)
        fastqc()                    Conduct FastQC analysis for barcode (sample)
        human_alignment()           Align reads to human reference genome using ont-specific parameters
        samtools_stats()            Perform samtools stats on the alignment output to calculate alignment stats
        hr_removal()                Convert sam to fastq (samtools), outputting only unmapped reads (non-human)
        classification()            Centrifuge classification. Create index summary file, run centrifuge
        create_kreport()            Convert centrifuge report to a kreport (kraken-style report)
        de_novo_assembly()          De novo assembly with metaFlye
        assembly_split()            Split assembly output into a file per assembly using faidx
        get_passing_contigs()       Parse assembly info file to get a list of contigs that pass minimum coverage
                                    requirements (60X)
        concatenate_sequences()     Concatenate input sequences into one fna file
        rename_headers()            Add run and barcode name into the read headers in fasta files
    """

    def __init__(self, fastq, runfolder, run_id, barcode, run_barcode, index_dir,
                 split_barcodes_dir, runfolder_out, logfile, reference_genome):
        """
        The constructor for FastqAnalysis class
        """
        # DATA INPUT PATHS
        self.run_id = run_id
        self.runfolder = runfolder
        self.barcode = barcode
        self.run_barcode = run_barcode
        self.fastq = fastq
        self.index_dir = index_dir
        self.index = "{}/p_compressed+h+v".format(self.index_dir)
        self.split_barcodes_dir = split_barcodes_dir
        self.reference_genome = reference_genome

        # DATA OUTPUT PATHS
        self.logfile = logfile
        # QC
        self.runfolder_out = runfolder_out
        self.qc_dir = "{}/QC".format(self.runfolder_out)
        self.fastqc_html_outfile = "{}/{}_fastqc.html".format(self.qc_dir, self.run_barcode)
        self.fastqc_zip_outfile = "{}/{}_fastqc.zip".format(self.qc_dir, self.run_barcode)
        # pycoqc
        self.sequencing_summary_file = "{}/sequencing_summary_{}.txt".format(self.split_barcodes_dir, self.barcode)
        self.pycoqc_outfile = "{}/{}_pycoQC_output.json".format(self.qc_dir, self.run_barcode)
        # centrifuge classification
        self.class_dir = "{}/centrifuge".format(self.runfolder_out)
        self.centrifuge_index_summary = "{}/centrifuge_index_summary.txt".format(self.runfolder_out)
        self.centrifuge_kreport = "{}/{}_kreport.tsv".format(self.class_dir, self.run_barcode)
        self.centrifuge_results = "{}/{}_results.tsv".format(self.class_dir, self.run_barcode)
        self.centrifuge_summary = "{}/{}_summary.tsv".format(self.class_dir, self.run_barcode)
        # human read removal
        self.hr_rem_dir = "{}/human_read_removal".format(self.runfolder_out)
        self.hr_rem_outpath = "{}/{}".format(self.hr_rem_dir, self.run_barcode)
        self.sam_output = "{}_aligned.sam".format(self.hr_rem_outpath)
        self.unmapped_fastq = "{}_unmapped.fastq".format(self.hr_rem_outpath)
        # samtools stats
        self.samtools_stats_file = "{}/{}_samtools_stats.txt".format(self.qc_dir, self.run_barcode)
        # metaflye de novo assembly
        self.assembly_dir = "{}/de_novo_assembly".format(self.runfolder_out)
        self.assembly_outdir = "{}/{}".format(self.assembly_dir, self.run_barcode)
        self.assembly_outfile = "{}/{}_assembly.fasta".format(self.assembly_outdir, self.run_barcode)
        self.assembly_log = "{}/flye.log".format(self.assembly_outdir)
        # Coverafge filtering
        self.assembly_split_dir = "{}/split_assemblies".format(self.runfolder_out)
        self.passing_contigs_dir = "{}/passing_contigs".format(self.runfolder_out)
        self.passing_contigs_fq = "{}/{}_passing_assembly.fasta".format(self.passing_contigs_dir, self.run_barcode)
        self.headers_renamed_fq = "{}_headers_renamed.fasta".format(self.passing_contigs_fq.split(".fasta")[0])


    def set_off_analysis(self):
        """
        Per-fastq analysis
        """
        # PRE-ALIGNMENT QC
        self.pycoqc()
        self.fastqc(self.fastq)

        # HUMAN READ REMOVAL
        # Remove human reads from the samples by alignment to the human reference genome (output sam file)
        self.human_alignment()
        # Samtools stats to show number of reads mapped to human genome
        self.samtools_stats()
        # Convert sam to fastq and remove human reads
        self.hr_removal()

        # FastQ on unmapped reads only (non-human reads only)
        self.fastqc(self.unmapped_fastq)

        # CLASSIFICATION
        self.classification()
        # convert classification output to centrifuge kreport
        self.create_kreport()

        # ASSEMBLY
        self.de_novo_assembly()

        # COVERAGE FILTERING
        # Split assembly output into a file per assembly using faidx
         # return a list of all filepaths containing passing contigs
        passing_contigs = self.assembly_split()
        # concatenate passing contigs into a single file
        self.concatenate_sequences(passing_contigs)
        # Add run and barcode name into the read headers in fasta files
        self.rename_headers()


    def pycoqc(self):
        """
        Create pycoQC json report for barcode (sample)
        """
        # If pycoQC not yet run, run pycoQC
        if os.path.isfile(self.sequencing_summary_file):
            if not os.path.isfile(self.pycoqc_outfile):
                logger(self.logfile).info("PYCOQC - Creating PycoQC json report for {}. ".format(self.run_barcode))
                logger(self.logfile).info("PYCOQC - version {}".format(app_dictionary["pycoqc"]))
                pycoqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} pycoQC -f {} --json_outfile " \
                                 "{}".format(app_dictionary["pycoqc"],
                                             self.sequencing_summary_file, self.pycoqc_outfile)
                run_command(self.logfile, pycoqc_command, "PYCOQC")
            else:
                logger(self.logfile).info("PYCOQC - Directory not empty - pycoQC "
                                          "already run for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("PYCOQC - Summary file for {} could not be found. ".format(self.run_barcode))


    def fastqc(self, fastq):
        """
        Conduct FastQC analysis for barcode (sample). Renames output files to ensure they match the naming conventions
        used in the analysis.

        NEEED TO CHECK THAT THSI RENAMES THE FILES CORRECTLY UNMAPPED AND NORMAL
        """
        outfile_pattern = fastq.split("/")[-1].split(".f")[0]
        old_html_outfile = "{}/{}_fastqc.html".format(self.qc_dir, outfile_pattern)
        old_zip_outfile = "{}/{}_fastqc.zip".format(self.qc_dir, outfile_pattern)

        # if output file already present, do not re-analyse
        if os.path.isfile(self.fastqc_html_outfile):
            logger(self.logfile).warning("FASTQC: Output for {} already exists".format(self.run_barcode))
        else:
            logger(self.logfile).info("FASTQC: Creating FastQC file for {}".format(self.run_barcode))
            logger(self.logfile).info("FASTQC - version {}".format(app_dictionary["fastqc"]))
            fastqc_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} " \
                             "fastqc {} -o {}".format(app_dictionary["fastqc"], fastq, self.qc_dir)
            run_command(self.logfile, fastqc_command, "FASTQC")

        # rename outfile to ensure it matches the naming conventions of outputs (so it appears in correct place in
        # MultiQC report)
        if os.path.isfile(self.fastqc_html_outfile):
            logger(self.logfile).warning("FASTQC: Html output with correct name for {} "
                                         "already exists".format(self.run_barcode))
        else:
            try:
                os.rename(old_html_outfile, self.fastqc_html_outfile)
            except:
                logger(self.logfile).error("FASTQC: FastQC html outfile {} could "
                                           "not be renamed".format(old_html_outfile))
            else:
                logger(self.logfile).info("FASTQC: FastQC html outfile {} successfully renamed".format(old_html_outfile))

        if os.path.isfile(self.fastqc_zip_outfile):
            logger(self.logfile).warning("FASTQC: zip output with correct name for {} "
                                         "already exists".format(self.run_barcode))
        else:
            try:
                os.rename(old_zip_outfile, self.fastqc_zip_outfile)
            except:
                logger(self.logfile).error("FASTQC: FastQC zip outfile {} could "
                                           "not be renamed".format(old_zip_outfile))
            else:
                logger(self.logfile).info("FASTQC: FastQC zip outfile {} successfully renamed".format(old_zip_outfile))


    def human_alignment(self):
        """
        Align reads to human reference genome using ont-specific parameters
        -a generates the outputs in the SAM format which allows extraction of unmapped reads
        -x -map-ont aligns long noisy reads of ~10% error rate to reference genome
        """

        if os.path.isfile(self.sam_output):
            logger(self.logfile).warning("MINIMAP2 HR MAPPING: Already complete for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("MINIMAP2 HR MAPPING: Aligning reads to human "
                                      "reference genome for {}".format(self.fastq))
            logger(self.logfile).info("MINIMAP2 - version {}".format(app_dictionary["minimap2"]))

            minimap2_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} minimap2 -ax " \
                               "map-ont {} {} -o {}".format(app_dictionary["minimap2"],
                                                            self.reference_genome, self.fastq, self.sam_output)
            run_command(self.logfile, minimap2_command, "MINIMAP2")


    def samtools_stats(self):
        """
        Perform samtools stats on the alignment output to calculate alignment stats
        (eg. no. reads mapped to human genome)
        """

        if os.path.isfile(self.samtools_stats_file):
            logger(self.logfile).warning("SAMTOOLS STATS: Samtools stats already "
                                         "conducted for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("SAMTOOLS STATS: Run samtools stats for {}".format(self.run_barcode))
            logger(self.logfile).info("SAMTOOLS - version {}".format(app_dictionary["samtools"]))

            samtools_stats_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} samtools stats " \
                                     "{} > {}".format(app_dictionary["samtools"],
                                                      self.sam_output, self.samtools_stats_file)
            run_command(self.logfile, samtools_stats_command, "SAMTOOLS STATS")


    def hr_removal(self):
        """
        Convert sam file to fastq file using samtools, outputting only unmapped reads (non-human)
        -f 4 = only outputs unmapped segments (with bitwise flag integer = 4)
        """
        if os.path.isfile(self.unmapped_fastq):
            logger(self.logfile).warning("SAMTOOLS FASTQ: Already complete for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("SAMTOOLS FASTQ: Convert non-human from sam "
                                      "to fastq for {}_unmapped.sam".format(self.hr_rem_outpath))
            logger(self.logfile).info("SAMTOOLS - version {}".format(app_dictionary["samtools"]))
            samtools_fastq_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} " \
                                     "samtools fastq -f 4 {} -0 {}".format(app_dictionary["samtools"],
                                                                            self.sam_output, self.unmapped_fastq)
            run_command(self.logfile, samtools_fastq_command, "SAMTOOLS FASTQ")


    def classification(self):
        """
        Centrifuge classification
        Create centrifuge index summary file, and run centriuge
        Disregard abundance - there is a bug in centrifuge that calculates the abundance as 0 for all taxIDs
        """
        # Create centrifuge index summary
        if os.path.isfile(self.centrifuge_index_summary):
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge index summary file already created")
        else:
            logger(self.logfile).info("CENTRIFUGE: Creating centrifuge index "
                                      "summary file for {}".format(self.run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(app_dictionary["centrifuge"]))
            # run centrifuge-inspect to output a summary of the index used for classification
            centrifuge_inspect = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} centrifuge-inspect -s " \
                                 "{} > {}".format(app_dictionary["centrifuge"], self.index,
                                                  self.centrifuge_index_summary)
            run_command(self.logfile, centrifuge_inspect, "CENTRIFUGE")

        # Run centrifuge classification
        if os.path.isfile(self.centrifuge_results) and os.path.isfile(self.centrifuge_summary):
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge already run for {}".format(self.run_barcode))
        else:
            # classify reads
            logger(self.logfile).info("CENTRIFUGE: Running centrifuge classification for {}".format(self.run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(app_dictionary["centrifuge"]))
            # -S is file to write classification results to, --report-file is file to write classification summary to,
            # -x is the index, -f is the input sequence, --env sets centriuge indexes environment variable so it knows
            # where to look for the index
            centrifuge_command = "docker run --env CENTRIFUGE_INDEXES={} --rm -v `pwd`:`pwd` -w `pwd` {} " \
                                 "centrifuge -x p_compressed+h+v -q {} -S {} " \
                                 "--report-file {}".format(self.index_dir, app_dictionary['centrifuge'],
                                                           self.unmapped_fastq, self.centrifuge_results,
                                                           self.centrifuge_summary)
            run_command(self.logfile, centrifuge_command, "CENTRIFUGE")


    def create_kreport(self):
        """
        Convert centrifuge report to a kreport (kraken-style report). This format is picked up by MultiQC, and can
        be loaded into Pavian
        """
        if os.path.isfile(self.centrifuge_kreport):
            logger(self.logfile).warning("CENTRIFUGE: Centrifuge-kreport already "
                                         "generated for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("CENTRIFUGE: Generating centrifuge-kreport for {}".format(self.run_barcode))
            logger(self.logfile).info("CENTRIFUGE - version {}".format(app_dictionary["centrifuge"]))
            kreport_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} centrifuge-kreport -x {}" \
                              " {} > {}".format(app_dictionary["centrifuge"], self.index,
                                                self.centrifuge_results, self.centrifuge_kreport)
            run_command(self.logfile, kreport_command, "CENTRIFUGE")


    def de_novo_assembly(self):
        """
        Runs using metaFlye.
            --meta enables mode for metagenome/uneven coverage assembly, designed for highly non-uniform coverage and
            is sensitive to underrepresented sequence at low coverage (as low as 2x)
            alternative haplotypes are collapsed
            minimum overlap length is chosen automatically based on read length distribution
            metaFlye performs one polishing iteration by default
        """
        # original non-sample named assembly output path that will be renamed
        assembly_outfile_old = "{}/assembly.fasta".format(self.assembly_outdir)

        if os.path.isfile(self.assembly_log):
            logger(self.logfile).warning("METAFLYE: De novo assembly already complete for {}".format(self.run_barcode))
        else:
            logger(self.logfile).info("METAFLYE: Performing METAFLYE de novo assembly for {}".format(self.unmapped_fastq))
            logger(self.logfile).info("METAFLYE - version {}".format(app_dictionary["flye"]))
            create_directory(self.assembly_outdir, self.logfile)
            metaflye_command = "docker run --rm -v `pwd`:`pwd` -w `pwd` {} flye --nano-raw {} --out-dir {} " \
                               "--meta".format(app_dictionary["flye"], self.unmapped_fastq, self.assembly_outdir)

            # rename output to add in run ID and barcode
            run_command(self.logfile, metaflye_command, "FLYE")

            if os.path.exists(assembly_outfile_old):
                try:
                    os.rename(assembly_outfile_old, self.assembly_outfile)
                except:
                    logger(self.logfile).error("METAFLYE: Assembly outfile {} could "
                                               "not be renamed".format(assembly_outfile_old))
                else:
                    logger(self.logfile).error("METAFLYE: Assembly outfile {} successfully "
                                               "renamed".format(assembly_outfile_old))
            else:
                logger(self.logfile).warning("METAFLYE: No assembly produced.")


    def assembly_split(self):
        """
        Split assembly output into a file per assembly using faidx.
        """
        # get dictionary of contigs that pass coverage requirements
        passing_contigs = self.get_passing_contigs(self.run_barcode, self.assembly_outdir)

        if os.listdir(self.assembly_split_dir):
            logger(self.logfile).warning("FAIDX: Assembly output already processed for {}".format(self.run_barcode))

        else:
            if os.path.exists(self.assembly_outfile):
                logger(self.logfile).info("FAIDX: Processing assembly output for {}".format(self.run_barcode))
                logger(self.logfile).info("PYFAIDX - version {}".format(app_dictionary["pyfaidx"]))
                faidx_command = "docker run -v `pwd`:`pwd` -v {}:{} {} /bin/sh -c 'cd {} ; " \
                                "faidx -x {}'".format(self.assembly_outdir, self.assembly_outdir,
                                                      app_dictionary["pyfaidx"], self.assembly_split_dir,
                                                      self.assembly_outfile)

                run_command(self.logfile, faidx_command, "FAIDX")
            else:
                logger(self.logfile).info("FAIDX: Assembly file not present for {}".format(self.run_barcode))

        # create dictionary to store only filepaths of contigs that surpass the coverage threshold of 60X
        assembly_list = []

        for file in os.listdir(self.assembly_split_dir):
            filepath = "{}/{}".format(self.assembly_split_dir, file)
            if any(contig in file for contig in passing_contigs):
                assembly_list.append(filepath)

        return assembly_list


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


    def concatenate_sequences(self, filepath_list):
        """
        Concatenate input sequences into one fna file
        """

        if filepath_list:
            if os.path.exists(self.passing_contigs_fq):
                logger(self.logfile).warning("CONCATENATE SEQUENCES: Passing sequencings already "
                                                  "concatenated for {}".format(self.run_barcode))
            else:
                logger(self.logfile).info("CONCATENATE SEQUENCES: Concatenating passing "
                                               "sequences for {}".format(self.run_barcode))
                assembly_string = ""
                for filepath in filepath_list:
                    assembly_string = "{} {}".format(assembly_string, filepath)

                concat_command = "cat {} > {}".format(assembly_string, self.passing_contigs_fq)
                run_command(self.logfile, concat_command, "CONCATENATE SEQUENCES")
        else:
            logger(self.logfile).warning("CONCATENATE SEQUENCES: No passing sequences "
                                         "found for {}".format(self.run_barcode))


    def rename_headers(self):
        """
        Add run and barcode name into the read headers in fasta files
        """

        if os.path.isfile(self.passing_contigs_fq):
            if os.path.isfile(self.headers_renamed_fq):
                logger(self.logfile).warning("RENAME HEADERS: Headers already renamed for {}".format(self.run_barcode))

            else:
                logger(self.logfile).warning("RENAME HEADERS: Renaming headers for {}".format(self.run_barcode))

                fasta = open(self.passing_contigs_fq, "r")
                newfasta = open(self.headers_renamed_fq, "w")

                for line in fasta.readlines():
                    if line.startswith(">"):
                        elements = line.split(">")
                        new_header = ">{}{}_{}".format(elements[0], self.run_barcode, elements[1])
                        newfasta.write(new_header)
                    else:
                        newfasta.write(line)
                fasta.close()
                newfasta.close()
        else:
            logger(self.logfile).warning("RENAME HEADERS: No passing sequences found for {}".format(self.run_barcode))


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

    logging.basicConfig(format='%(asctime)s %(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG,
                        handlers=[
                            logging.FileHandler(logfile),
                            logging.StreamHandler(sys.stdout)]
                        )
    return logging

def create_directory(directory, logfile):
    """
    Create output directories
    """
    if os.path.exists(directory):
        logger(logfile).warning("Directory not created - already exists: {}".format(directory))
    else:
        os.mkdir(directory)
        logger(logfile).info("Output directory created: {}".format(directory))
    return directory


def run_command(logfile, command, tool):
    """"
    Run the command as a subprocess using Popen, write the stdout and stderr to file using logging, and write to
    log dependent on exitcode
    """
    logger(logfile).info("{} - RUNNING COMMAND: {}".format(tool, command))
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
            logger(logfile).info('got line from subprocess: {}'.format(line))


def main():
    args = arg_parse()
    runfolder_dir = args['runfolder_dir']
    reference_dir = args['ref_dir']
    index_dir = args['centrifuge_index']
    out_dir = args['out_dir']
    ref_gen = args['ref_genome']

    all_runs_analysis = AllRuns(runfolder_dir, reference_dir, index_dir, out_dir, ref_gen)
    all_runs_analysis.set_off_analysis()


if __name__ == '__main__':
    main()