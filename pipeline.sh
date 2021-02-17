#!/bin/bash
# Install containers
python3 install_containers.py
# Load singularity
module load apps/singularity

# Create output directories
output_directories=(
    fastqc \
    pycoqc \
    human_read_removal
)

for i in "${output_directories[@]}"; do
    if [ -d "output/${i}" ];then
        echo "${i} directory already exists"
    else
        mkdir output/${i}
        echo "output folder ${i} created"  
    fi
done

# DATA ANALYSIS ==============================================================================================
# FASTQC -----------------------------------------------------------------------------------------------------
# for each fastq file
for i in S54 S62; do
  file="${i}_24hrs_fastqc.html"
  # if output file already present, do not re-analyse
  if [ -f output/fastqc/${file} ] ; then
    echo "FastQC output for ${file} already exists"
  else
    echo "Creating fastqc file for ${file}"
    singularity exec apps/fastqc.sif fastqc data/${i}_24hrs.fq -o output/fastqc
  fi
done

# PYCOQC ------------------------------------------------------------------------------------------------------
#PycoQC, with guppy barcoding file
# if not yet split, split the summary sequencing files according to barcodes
if [ ! "$(ls -A output/pycoqc)" ] ; then
  echo "Summary sequencing files not yet split"
  for file in data/sequencing_summary_*; do
    echo "Splitting summary sequencing file ${file} according to barcodes"
    singularity exec apps/pycoqc.sif Barcode_split --output_unclassified --min_barcode_percent 0.0 --summary_file ${file} --output_dir output/pycoqc/
  done
else
    echo "Directory not empty - barcodes already split"
fi
# Create pycoQC json report per barcode
for file in output/pycoqc/sequencing_summary_*; do
  #get barcode name
  barcode=$(echo "$file" | cut -d '_' -f 3 | cut -d '.' -f 1)
  echo "Creating pycoQC json report for ${file}"
  singularity exec apps/pycoqc.sif pycoQC -f ${file}  --json_outfile output/pycoqc/${barcode}_pycoQC_output.json
  done

# HUMAN READ REMOVAL -------------------------------------------------------------------------------------------
ref=data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna
for file in data/*.fq; do
  # select run  name
  run_name=$(echo "$file" | cut -d '.' -f 1 | rev | cut -d '/' -f 1 | rev)
  # align reads to human reference genome
  echo "Aligning reads to human reference genome for ${file}"
  singularity exec apps/minimap2.sif minimap2 -ax map-ont ${ref} ${file} > output/human_read_removal/${run_name}_aligned.sam
  # export unassigned reads to bam file with samtools
  #singularity shell apps/samtools -c "samtools view -f 4 file.bam > unmapped.sam"
done
# align reads using ont-specific parameters

# minimap2 for mapping alignment, bcftools consensus generation, SNP-sites to identify SNPs between samples, multi-locus sequence typing using MLST, and SNP-dists to calculate SNP distances.1

# MULTIQC ------------------------------------------------------------------------------------------------------
# create multiqc report, pulling in outputs from other tools
echo "Creating MultiQC report for the analysis"
singularity exec apps/multiqc.sif python -m multiqc output --outdir output/multiqc
