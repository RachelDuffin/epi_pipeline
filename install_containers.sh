# load singularity on HPC
# module load apps/singularity
# download containers (e.g. medaka)

cd apps

echo "INSTALLING QC TOOLS"
echo "Installing FastQC"
singularity pull --name fastqc docker://quay.io/biocontainers/fastqc:0.11.9--0
echo "Installing MultiQC"
singularity pull --name multiqc docker://quay.io/biocontainers/multiqc:1.9--py_1
echo "Installing PycoQC"
singularity pull --name pycoqc docker://quay.io/biocontainers/pycoqc:2.5.2--py_0


echo "INSTALLING PROCESSING TOOLS"
echo "Installing Minimap2"
singularity pull --name minimap2 docker://quay.io/biocontainers/minimap2:2.17--hed695b0_3
echo "Installing Samtools"
singularity pull --name samtools docker://quay.io/biocontainers/samtools:0.1.18--hfb9b9cc_10
echo "Installing bam2fastx"
singularity pull --name bam2fastx docker://quay.io/biocontainers/bam2fastx:1.3.1--he1c1bb9_0
echo "Installing Medaka"
singularity pull --name medaka docker://quay.io/biocontainers/medaka:1.2.1--py38hfcf0ad1_0

cd ..