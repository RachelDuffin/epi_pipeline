# Activate virtual environment
source venv/bin/activate
# Install singularity
./install_singularity.sh
# Install containers
./install_containers.sh

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
    singularity shell apps/fastqc -c "fastqc data/${i}_24hrs.fq -o output/fastqc"
  fi
done

# PYCOQC ------------------------------------------------------------------------------------------------------
#PycoQC, with guppy barcoding file
# split summary sequencing files according to barcodes
if [ ! "$(ls -A output/pycoqc)" ] ; then
  singularity shell apps/pycoqc -c "Barcode_split --output_unclassified --min_barcode_percent 0.0 --summary_file data/sequencing_summary_FAO06374_2bff58da.txt --output_dir output/pycoqc/"
for file in output/pycoqc/sequencing_summary_*; do
  echo $file
  #get barcode name
  barcode=$(echo $file | cut -d '_' -f 3 | cut -d '.' -f 1)
  #create pycoQC html report per barcode
  singularity shell apps/pycoqc -c "pycoQC -f ${file}  --json_outfile output/pycoqc/${barcode}_pycoQC_output.json"
  done
else
  echo "Directory not empty - barcodes already split"
fi

singularity shell apps/multiqc -c "python -m multiqc output --outdir output/multiqc"
pkill singularity
# create pycoqc report
#singularity shell apps/pycoqc -c "pycoQC -b data/sequencing_summary_FAO06374_2bff58da.txt " \
 #                                 "–f data/final_summary_FAO06374_2bff58da.txt " \
  #                                "-o output/pycoQC_output.html -j output/pycoQC_output.json"