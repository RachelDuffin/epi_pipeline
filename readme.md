# Outbreak Analysis Pipeline

This pipeline (pipeline.py) was created as part of an MSc project. It was written to be run on unix-based systems.  

The pipeline contains the following steps:

| Step                                       | Function                                                                                                                                                                                           | Tool                                               |
|--------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------|
| Pre-alignment QC                           | Generation of quality control metrics and plots from sequence data pre and post processing                                                                                                         | FastQC v0.11.9<br/>pycoQC v2.5.2                   |
 | Mapping alignment to human genome          | Mapping reads to the host genome (GRCh38) to allow human read removal                                                                                                                              | minimap2 v2.17-r941                                |
| Post-alignment QC                          | Collects alignment statistics from BAM files                                                                                                                                                       | samtools stats (samtools v1.12 Using htslib v1.12) |
| Human read removal                         | Computationally subtracting human reads from the dataset                                                                                                                                           | samtools fastq (samtools v1.12 Using htslib v1.12) |
| Classification                             | Provides a record of all species present within the sample                                                                                                                                         | centrifuge v1.0.4                                  |
| K-report generation                        | Generate a kreport which is parsable by MultiQC and can be viewed in Pavian                                                                                                                        | centrifuge v1.0.4                                  |
| De novo assembly (24 hour analysis)        | Assembling reads with no prior knowledge of the correct order of fragments, after 24 hours of sequencing. The assembly is required to conduct subsequent stages of the workflow including typing.  | metaFlye v2.8.3-b1695                              |
 | Filter assemblies by coverage              | Remove contigs <60X, as for de novo assembly a raw read coverage of 50-60X is required to generate enough reads to cover repetitive regions in the assembly                                        | pyfaidx v0.6.4, and self-written python code       |
| Add run and barcode name into contig names | Allows assemblies to be traced back to run, barcode, and so patient                                                                                                                                | Self-written python code                           |
| Multi locus sequence typing (MLST)         | Use the DNA sequences of internal fragments of multiple housekeeping genes to characterise isolates of microbial species.                                                                          | mlst v2.19.0                                       |
 | Aggregate QC                               | Aggregates results for all samples in analysis, allowing comparison and scrutinisation of quality metrics and results.                                                                             | MultiQC v1.12                                      |

Examples of some of the outputs from two analyses created by running this pipeline can be found in the pipeline_outputs 
directory, for a set of metagenomic samples and a set of monomicrobial samples. This includes the Krona reports 
generated in downstream analyses. Just follow these links to view:
* [Metagenomic MultiQC Report](https://rachelduffin.github.io/epi_pipeline/pipeline_outputs/clinical_metagenomic_samples/metagenomic_multiqc_report.html)
* [Monomicrobial MultiQC Report](https://rachelduffin.github.io/epi_pipeline/pipeline_outputs/monomicrobial_samples/monomicrobial_multiqc_report.html)
* [Metagenomic Krona Report](https://rachelduffin.github.io/epi_pipeline/pipeline_outputs/clinical_metagenomic_samples/metagenomic_taxonomy_krona.html)
* [Monomicrobial Krona Report](https://rachelduffin.github.io/epi_pipeline/pipeline_outputs/monomicrobial_samples/monomicrobial_taxonomy_krona.html)

All tools used in the pipelines are biocontainers pulled and run using Docker, with specific SHA tags for each viewable 
in the config file (config.py)

## Setup 
Once this repository has been pulled, navigate into the project directory, and:
1. Create a virtual environment using:
```
pip3 install virtualenv
virtualenv venv 
source venv/bin/activate
```
2. Install the requirements.txt file:
```
pip3 -r install requirements.txt
```
## Required inputs, and running the pipeline
The pipeline.py script must be run from within the same directory as the script, in and all scripts and 
files should be left in their original positions within the directory.

The script takes several command line arguments:
```
-d      --runfolder_dir         Directory containing runfolders to include in the analysis
-r      --ref_dir               Directory containing reference sequences
-i      --centrifuge_index      Directory containing centrifuge index
-o      --out_dir               Directory to create output directories within
-g      --ref_genome            Human reference genome fna file
```
For example, it could be run as follows:
```
python3 pipeline.py -d /home/dir_containing_runfolders/ -r /home/reference_sequences -i /home/centrifuge_index/ -o /home/output_dir/ -g /home/GCF_000001405.39_GRCh38.p13_genomic.fna
```

Input run folders should have the following directory structure (i.e. a run folder per run within the run folders 
directory, containing the final summary and sequencing summary files, as well as a .fq file per barcode).
```
dir_containing_runfolders/
    runfolder_1/
        barcode_1.fq
        barcode_2.fq
        barcode_3.fq
        etc.
        sequencing_summary_*.txt
        final_summary_*.txt 
    /runfolder_2/
        etc.
```
Once within the directory containing the "pipeline.py" script, the following code should be run:

```
python3 pipeline.py
```

## Outputs 
The script will create a directory per run within the output directory. Within each of these directories, 
a subdirectory containing the output from each tool will be created. 

Resulting output files can then be copied off the cluster using the following command on your personal linux computer:

```
rsync -avz -partial USERNAME@login.rosalind.kcl.ac.uk:REMOTE_DIRECTORY_PATH LOCAL_DIRECTORY_PATH
```

The centrifuge k-reports for all samples in the analysis can be loaded into a Pavian webserver which can be spun up 
using docker with the following commands: 
```
docker pull 'florianbw/pavian'
docker run --rm -p 5000:80 florianbw/pavian
```
This webserver can then be accessed at 0.0.0.0:5000 and the kreports loaded in for interactive visualisation of 
taxonomic abundance within each sample, in the form of tables and per-sample sankey visualisations.

Krona can be used to explore relative taxonomic abundance within each sample, in the form of radial space-filling 
hierarchical visualizations. This can be achieved by running your modification of the below commands, the first on the 
command line and the second within the interactive container:
```
sudo docker run -it --rm -v /home/rachel:/home/rachel quay.io/biocontainers/krona@sha256:8917b9840b369d102ee759a037cc8577295875952013aaa18897c00569c9fe47 /bin/bash 
apt-get update && apt-get install make && ktUpdateTaxonomy.sh && cd /home/rachel && ktImportTaxonomy /home/rachel/msc_project/cut_summary_files/metagenomic/*_cut.tsv
```
This will generate an interactive krona.html file which can be opened in a web browser.

The multiqc.html report can be opened in web browser to view the outputs. 

## Assembly testing scripts
The test_scripts and assembly_results directories contain data relating to the assembler tool tests that were run as 
part of the project. 
* assembly_results contains csv files containing the results of performance and accuracy tests, the R script written to
create plots from these results, and .png files of the plots created
* test_scripts contains the scripts created to conduct the performance and accuracy tests (these were run on the 
Rosalind High Performance Cluster):

| Script                        | Function                                                                                                                                                                                                                          | 
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| test_assemblers.py            | Testing performance of canu, raven and flye (hard-coded filepaths)                                                                                                                                                                |
| test_accuracy.py              | Compare assembly outputs with reference sequences to determine accuracy (command line input filepaths)                                                                                                                            |
| generate_simulated_data.py    | Incomplete, written to generate simulated monomicrobial and metagenomic datasets from reference genomes for testing datasets using NanoSim                                                                                        |
| count_reads.py                | Calculates number of reads and lengths of reads in fastq files and outputs the metrics to text files. This was written to calculate metrics to input into NanoSim commands. Takes input filepaths from a dictionary (hard-coded). |

The summary image of all the results generated using the test_assemblers script for performance testing of Canu, 
Raven and metaFlye is displayed below:

![alt text](https://github.com/RachelDuffin/epi_pipeline/blob/master/assembly_results/Assembly_test_plots.png?raw=true)