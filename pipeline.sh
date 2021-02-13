# Activate virtual environment
# Install singularity
sudo ./install_singularity.sh
# Install containers
./install_containers.sh

source venv/bin/activate
sudo singularity config global --set 'mount tmp' no

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

# get current working directory
wd=$(pwd)

#set singularity bind paths to bind data on the host to the containers
#export SINGULARITY_BIND="${wd}/output,${wd}/data,${wd},${wd}/output/fastqc"
export SINGULARITY_SHELL=/bin/bash
export LANG=en_GB.UTF-8

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
    sudo singularity exec apps/fastqc.sif fastqc data/${i}_24hrs.fq -o output/fastqc
  fi
done

# PYCOQC ------------------------------------------------------------------------------------------------------
#PycoQC, with guppy barcoding file
# split summary sequencing files according to barcodes
if [ ! "$(ls -A output/pycoqc)" ] ; then
  singularity exec apps/pycoqc Barcode_split --output_unclassified --min_barcode_percent 0.0 --summary_file data/sequencing_summary_FAO06374_2bff58da.txt --output_dir output/pycoqc/ --force
for file in output/pycoqc/sequencing_summary_*; do
  #get barcode name
  barcode=$(echo "$file" | cut -d '_' -f 3 | cut -d '.' -f 1)
  #create pycoQC json report per barcode
  singularity exec apps/pycoqc.sif pycoQC -f ${file}  --json_outfile output/pycoqc/${barcode}_pycoQC_output.json
  done
else
  echo "Directory not empty - barcodes already split"
fi

# HUMAN READ REMOVAL -------------------------------------------------------------------------------------------
ref=/home/rachel/outbreak_pipeline/data/human_genome/ncbi/GCF_000001405.39_GRCh38.p13_genomic.fna
for file in /home/rachel/outbreak_pipeline/data/sample_fasta/*; do
  # select run  name
  run_name=$(echo "$file" | cut -d '.' -f 1 | rev | cut -d '/' -f 1 | rev)
  # align reads to human reference genome
  singularity exec apps/minimap2.sif minimap2 -ax map-ont $ref $file > output/human_read_removal/${run_name}_aligned.sam
  # export unassigned reads to bam file with samtools
  #singularity shell apps/samtools -c "samtools view -f 4 file.bam > unmapped.sam"
done
# align reads using ont-specific parameters


# MULTIQC ------------------------------------------------------------------------------------------------------
# create multiqc report, pulling in outputs from other tools
singularity exec apps/multiqc.sif python -m multiqc output --outdir output/multiqc
