# Outbreak Analysis Pipeline

This pipeline (pipeline.py) was created as part of an MSc project. It contains the following steps:

| Step | Function | Tool |
| ----------|-----------|-----------|
| QC | Generation of quality control metrics and plots from sequence data pre and post processing | |
| Human read removal | Mapping reads to the host genome (GRCh38) and computationally subtracting these from the dataset | |
| De novo assembly (24 hour analysis) | Assembling reads with no prior knowledge of the correct order of fragments, after 24 hours of sequencing | |
| Mapping assembly (6 hour analysis)| | |

## Requried inputs, and running the pipeline

The pipeline has been written to be run on the Rosalind KCL high performance cluster. Once the code has been pulled, the following directory structure should be replicated:

```
/outbreak_pipeline
    pipeline.py
    install_containers.py
    readme.md
    /data
        /human_genome
        /run_folders
            */run_name_1
                final_summary_{}.txt 
                barcode_1.fq
                barcode_2.fq
                barcode_3.fq
            /run_name_2
                etc.
    /output
```
This means creating a run folder per run within the run_folders directory, which should contain the final_summary text file, and a fastq file per barcode from that run (containing the concatenated reads).

Once within the directory
containing the "pipeline.py" script, the following code should be run:

```
python3 pipeline.py
```

## Outputs 
The script will create a directory per run within the outputs directory. Within each of these directories, a subdirectory containing the output from each tool will be created. Among these outputs is a MultiQC file per run. 

Resulting output files can then be copied off the cluster using the following command on your personal linux computer:

```
rsync -avz -partial USERNAME@login.rosalind.kcl.ac.uk:REMOTE_DIRECTORY_PATH LOCAL_DIRECTORY_PATH
```

## Assembly testing scripts
The scripts that were created for testing tools are as follows:
* test_assemblers.py - tests performance of canu, raven and flye (hard-coded filepaths)
* test_accuracy.py - compares assembly outputs with reference sequences to determine accuracy (command line 
input filepaths)
* count_reads.py - Calculates number of reads and lengths of reads in fastq files and outputs the metrics to text files. 
Takes input filepaths from a dictionary (hard-coded). 
* generate_simulated_data.py - generates simulated monomicrobial and metagenomic datasets for testing purposes from 
reference genomes (in a reference sequence directory). 
* 
