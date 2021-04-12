import os
import subprocess
# pull singularity images (e.g. medaka)

app_dictionary = {
    "fastqc": "quay.io/biocontainers/fastqc@sha256:319b8d4eca0fc0367d192941f221f7fcd29a6b96996c63cbf8931dbb66e53348",
    "multiqc": "quay.io/biocontainers/multiqc@sha256:2c3362c12409ca8f06f182c35cf8a804ee23cf8a32dc02c952fa60bc6e48db96",
    "pycoqc": "quay.io/biocontainers/pycoqc@sha256:ea0a084751a0b48b5ffe90e9d3adfa8f57473709a1b0a95c9cb38d434ee3a9a2",
    "minimap2":
        "quay.io/biocontainers/minimap2@sha256:7f95eecc8eeee8ef8ae7e24d1d1a49ee056fb21d72aea4e2def97484f8a206c5",
    "samtools":
        "quay.io/biocontainers/samtools@sha256:2b911396b907769945a5446c8779d7be83f999bb4211733a8259040e67f4065f",
    "bam2fastx":
        "quay.io/biocontainers/bam2fastx@sha256:c67f2281e6995edf2525e22271f8b68c483b2bb2866c5b66f7eec519348e5312",
    "medaka": "quay.io/biocontainers/medaka@sha256:6aa52d718af0f48cf6630e88b22cd7187770bbf028ef89bc54ec7fad2ff7a35f",
    "flye": "quay.io/biocontainers/flye@sha256:f895c72298ea3ae568c265cfb575fefeca768c42870cfea0ef3a4cfc35704086",
    "mlst": "quay.io/biocontainers/mlst@sha256:d8d8e731f165df398d95ad0969333aef45bcae678da11a26a1d2e76d44bc2698",
    "snp-sites":
        "quay.io/biocontainers/snp-sites@sha256:02c868581e36f97cc16f066f0ef4c2d169873ca4cf6932a104beb10c828a9c5c",
    "snp-dists":
        "quay.io/biocontainers/snp-dists@sha256:7ac5037a7967252593acee51b2ce2202b30a387d53bd62005fc3e927145557b4",
    "pyfaidx": "quay.io/biocontainers/pyfaidx@sha256:e7a2a24d3b5ea10181be085e0b39056c8c90a3298ab13f3f33104d2e3f46f314"
}


def install_tools():
    """
    Checks if each docker image specified in app_dictionary already exists, if not pulls from biocontainers
    """
    print("--------------------------\nInstalling pipeline tools\n--------------------------")

    for key in app_dictionary:
        if subprocess.check_output("sudo docker images -q " + app_dictionary[key], shell=True):
            print(key + " docker image already pulled")
        else:
            print("Installing " + key)
            os.system("sudo docker image pull " + app_dictionary[key])
    print("----------------------------")


def main():
    install_tools()


if __name__ == '__main__':
    main()
