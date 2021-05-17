import install_containers
import os
import subprocess

app_dictionary = {"centrifuge":
        "quay.io/biocontainers/centrifuge@sha256:e3ce6d3d83a1df5327ee27b66d4c9eedb05b7bcd2eae8f78a7f7b9c1e8672c1c"
                  }

def main():
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])
    download_taxonomy = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} centrifuge-download -o taxonomy " \
                        "taxonomy".format(app_dictionary["centrifuge"])
    download_genomes = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} centrifuge-download -o library -m -d " \
                        "'archaea,bacteria,viral' refseq > seqid2taxid.map".format(app_dictionary["centrifuge"])
    concatenate_sequences = "cat library/*/*.fna > input-sequences.fna"

    build_index = "sudo docker run --rm -v `pwd`:`pwd` -w `pwd` -i -t {} centrifuge-build -p 4 --conversion-table " \
                  "seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp " \
                  "input-sequences.fna abv".format(app_dictionary["centrifuge"])

    for docker_command in download_taxonomy, download_genomes, concatenate_sequences, build_index:
        subprocess.run(docker_command, shell=True)


if __name__ == '__main__':
    main()
