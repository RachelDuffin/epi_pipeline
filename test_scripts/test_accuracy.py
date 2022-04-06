"""
test_accuracy: Script to compare accuracy between assembly tools. N.B. this script will only work in the Rosalind
               HPC environment.

Takes command line arguments as inputs:
    -r = Input reference sequence file against which to compare your query sequence
    -q = Input query sequence to compare to the reference
    -o = Directory to save file outputs

ACCURACY:
Comparing the assembly to the reference genome, using the MUMmer package. Compare by:
    - No. bases in assembly (should be length of ref. genome)
    - No. contigs (should be 1)
    - Average identity - % similarity between assembly and reference genome (should be near to 100%)
    - Coverage - ratio of no. aligned bases in ref to length of ref (should be near to 100%)
    - No. mismatches - no. single-base differences between assembly and ref (should be 0)
    - No. indels - number of insertions and deletions between assembly and ref

"""
import argparse
import os
import subprocess
import shlex
from pathlib import Path
import test_assemblers
import glob


tool_dict = {
    "mummer" : "quay.io/biocontainers/mummer4@sha256:71a50b71457f7092177508522cd48da25ab0ece07544b946756b47f0dbfbff39"
}

pattern_dict = {
    "flye": "*assembly.fasta",
    "raven": "*stdout_*_raven_*_thread.txt",
    "canu": "*.contigs.fasta"
}

def arg_parse():
    """
    Parses arguments supplied by the command line.
        :return: (Namespace object) parsed command line attributes
    Creates argument parser, defines command line arguments, then parses supplied command line arguments using the
    created argument parser.
    """
    parser = argparse.ArgumentParser(description="Read file paths from command line")
    parser.add_argument('-r', '--ref_sequence', dest = "ref_filename", required = True,
                        help = "Input reference sequence file against which to compare your query sequence",
                        metavar = "FILE")
    parser.add_argument('-q', '--query_sequence', dest= "query_filename", required = True,
                        help = "Input query sequence to compare to the reference")
    parser.add_argument('-o', '--out_dir', dest= "out_dir", required = True,
                        help = "Directory to save file outputs")
    return parser.parse_args()


def run_comparison(reference_file, query_file, out_dir, assembler, base_path):
    """
    Compare accuracy of the output fasta files to the reference genome using dnadiff:
        - Number of bases (should equal length of reference)
        - Number of contigs (less the better)
        - Average identity (higher is better)
        - Coverage (higher is better)
        - Number of mismatches (less the better)
        - Number of indels (less the better)
        OUTPUT:
    """
    cwd = os.getcwd()
    os.chdir(out_dir)

    print("RUNNING {} COMPARISON: FOR {} AND {}".format(assembler, reference_file, query_file))

    command = 'module load apps/singularity && singularity exec --bind `pwd`:`pwd` ' \
              '{}/mummer.sif dnadiff {} {}'.format(base_path, reference_file, query_file)

    # run comparison
    process = subprocess.Popen(command, shell=True, universal_newlines=True, bufsize=1,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # collect stdout and stderr and print
    stdout, stderr = process.communicate()
    print(stdout, stderr)

    print("EXECUTED COMMAND: " + command)

    os.chdir(cwd)

def get_file_list(path, assembler, pattern):
    """
    Get list of files matching pattern
    """
    os.chdir(path)
    file_list = glob.glob("{}/*/{}/*/{}".format(path, assembler, pattern), recursive=True)
    return file_list

def main():
    args = arg_parse()
    print(args)
    out_dir = args.out_dir
    ref = args.ref_filename
    query = args.query_filename

    # path the singularity image is saved to
    base_path = os.getcwd().rsplit('/', 2)[0]

    # get singularity image
    for tool in tool_dict:
        test_assemblers.singularity_pull(tool, tool_dict[tool])

    # run comparisons for each assembler tool
    assemblers = ["flye", "raven", "canu"]

    for assembler in assemblers:
        # create output accuracy directory if doesn't already exist
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        # run comparison for input query sequence
        run_comparison(ref, query, out_dir, assembler, base_path)

if __name__ == '__main__':
    main()

