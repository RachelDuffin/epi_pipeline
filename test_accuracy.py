"""
test_accuracy: Script to compare accuracy between assembly tools

ACCURACY:
Comparing the assembly to the reference genome, using the MUMner package. Compare by:
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

app_dictionary = {
    "mummer" : "quay.io/biocontainers/mummer4@sha256:71a50b71457f7092177508522cd48da25ab0ece07544b946756b47f0dbfbff39"
}

def arg_parse():
    """
    Parses arguments supplied by the command line.
        :return: (Namespace object) parsed command line attributes
    Creates argument parser, defines command line arguments, then parses supplied command line arguments using the
    created argument parser.
    """
    parser = argparse.ArgumentParser(description="Read file paths from command line")
    parser.add_argument('-r', '--ref_sequence', dest = "ref_filename", required = True, type = validate_file,
                        help = "Input reference sequence file against which to compare your query sequence",
                        metavar = "FILE")
    parser.add_argument('-q', '--query_sequence', dest= "query_filename", required = True, type = validate_file,
                        help = "Input query sequence to compare to the reference")
    parser.add_argument('-o', '--out_dir', dest= "out_dir", required = True, type = validate_file,
                        help = "Directory to save file outputs")
    return parser.parse_args()


def run_comparison(reference_file, query_file, out_dir):
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
    chdir(out_dir)
    command = 'dnadiff {} {}'.format(reference_file, query_file)

    print("Executed command: " + command)
    subprocess.check_output(shlex.split(command))
    chdir(cwd)


def main():
    args = arg_parse()
    run_comparison(args.ref_filename, args.query_filename, args.out_dir)

if __name__ == '__main__':
    main()

