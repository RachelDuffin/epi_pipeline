import os
# load singularity on HPC
# module load apps/singularity
# download containers (e.g. medaka)

app_dictionary = {
    "fastqc.sif": "docker://quay.io/biocontainers/fastqc:0.11.9--0",
    "multiqc.sif": "docker://quay.io/biocontainers/multiqc:1.9--py_1",
    "pycoqc.sif": "docker://quay.io/biocontainers/pycoqc:2.5.2--py_0",
    "minimap2.sif": "docker://quay.io/biocontainers/minimap2:2.17--hed695b0_3",
    "samtools.sif": "docker://quay.io/biocontainers/samtools:0.1.18--hfb9b9cc_10",
    "bam2fastx.sif": "docker://quay.io/biocontainers/bam2fastx:1.3.1--he1c1bb9_0",
    "medaka.sif": "docker://quay.io/biocontainers/medaka:1.2.1--py38hfcf0ad1_0",
    "flye.sif" : "docker://quay.io/biocontainers/flye:2.8.3--py36h5202f60_0"
}
save_location = "apps"


def install_tools():
    print("--------------------------\nInstalling pipeline tools")
    for key in app_dictionary:
        if os.path.exists(save_location + "/" + key):
            print(key + " already installed")
        else:
            print("Installing " + key)
            os.system("module load apps/singularity; singularity pull --dir " + save_location + " " + key + " " +
                      app_dictionary[key])
    print("----------------------------")
