"""
Installs containers required in the pipeline
"""

import os
import subprocess

app_dictionary = {
    "fastqc": "quay.io/biocontainers/fastqc@sha256:319b8d4eca0fc0367d192941f221f7fcd29a6b96996c63cbf8931dbb66e53348",
    "multiqc": "quay.io/biocontainers/multiqc@sha256:2c3362c12409ca8f06f182c35cf8a804ee23cf8a32dc02c952fa60bc6e48db96",
    "pycoqc": "quay.io/biocontainers/pycoqc@sha256:ea0a084751a0b48b5ffe90e9d3adfa8f57473709a1b0a95c9cb38d434ee3a9a2",
    "minimap2":
        "quay.io/biocontainers/minimap2@sha256:7f95eecc8eeee8ef8ae7e24d1d1a49ee056fb21d72aea4e2def97484f8a206c5",
    "samtools":
        "quay.io/biocontainers/samtools@sha256:2b911396b907769945a5446c8779d7be83f999bb4211733a8259040e67f4065f",
    "flye": "quay.io/biocontainers/flye@sha256:f895c72298ea3ae568c265cfb575fefeca768c42870cfea0ef3a4cfc35704086",
    "mlst": "quay.io/biocontainers/mlst@sha256:d8d8e731f165df398d95ad0969333aef45bcae678da11a26a1d2e76d44bc2698",
    "mlstcheck" : 
        "quay.io/biocontainers/perl-bio-mlst-check@sha256:"
        "828f4a536603a39559118a279e9c66d71e83a155cf0ac824efdb9339ba59e201",
    "snp-sites":
        "quay.io/biocontainers/snp-sites@sha256:02c868581e36f97cc16f066f0ef4c2d169873ca4cf6932a104beb10c828a9c5c",
    "snp-dists":
        "quay.io/biocontainers/snp-dists@sha256:7ac5037a7967252593acee51b2ce2202b30a387d53bd62005fc3e927145557b4",
    "centrifuge":
        "quay.io/biocontainers/centrifuge@sha256:c43220b8fc171eb4a1b88fc1d6df104c51172a8463bb96c10b03cb315f41fb2e"
}


def install_tools(key, value):
    """
    Checks if docker image specified in app_dictionary already exists, if not pulls from biocontainers
    """
    if subprocess.check_output("docker images -q " + value, shell=True):
        print("INSTALLING DOCKER IMAGE: " + key + " docker image already pulled")
    else:
        print("INSTALLING DOCKER IMAGE: " + "Installing " + key)
        os.system("docker image pull " + value)
    print("----------------------------")


def main():
    for key in app_dictionary:
        value = app_dictionary[key]
        install_tools(key, value)


if __name__ == '__main__':
    main()
