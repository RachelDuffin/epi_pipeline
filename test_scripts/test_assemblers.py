"""
Script to compare performance, scalability and accuracy between assembly tools.

Conduct comparison of metagenomic assemblers, using 5 E. faecium isolates (monomicrobial samples), and Zymo mock
community data

PERFORMANCE AND SCALABILITY:
Comparison of performance and memory usage between tools (single thread), using /usr/bin/time -v. Compare by:
- Wall clock time - elapsed time from start to end of program
- CPU time - Time the CPU spends in user mode (running program's code) and kernel mode (executing system calls)
- Peak memory usage - max memory used by program during lifetime
- Parallel speedup - ratio of time to run program with one thread to time to run with N threads
This is run multiple times changing the number of threads and observing corresponding changes in metrics
(all tools - flye, canu, raven - have multithreading support). The linux box has 20 cores.
"""
import subprocess
import os
import sys

assembler_dictionary = {
    "raven":
        "quay.io/biocontainers/raven-assembler@sha256:3bc4cc61483cc48243f6b416eaae41f24eb95f75b7a2770f8062c75b5ac53da3",
    "flye": "quay.io/biocontainers/flye@sha256:f895c72298ea3ae568c265cfb575fefeca768c42870cfea0ef3a4cfc35704086",
    "canu": "quay.io/biocontainers/canu@sha256:b48b52afc355477015ef60bebded3b4ab3d3099bbf5698879de8eb600c9ff1a4"
}

input_filepaths = {
    "/data/monomicrobial_samples/enterococcus_faecium/ef1_bc_75/210612_EF_R1_barcode01.fq":
        "/output/assemblers/monomicrobial_samples/",
    "/data/monomicrobial_samples/enterococcus_faecium/ef1_bc_75/210612_EF_R1_barcode02.fq":
        "/output/assemblers/monomicrobial_samples/",
    "/data/monomicrobial_samples/enterococcus_faecium/ef1_bc_75/210612_EF_R1_barcode03.fq":
        "/data/monomicrobial_samples/assemblers/monomicrobial_samples/",
    "/data/monomicrobial_samples/enterococcus_faecium/ef1_bc_75/210612_EF_R1_barcode04.fq":
        "/output/assemblers/monomicrobial_samples/",
    "/data/monomicrobial_samples/enterococcus_faecium/ef1_bc_75/210612_EF_R1_barcode05.fq":
        "/output/assemblers/monomicrobial_samples/",
    "/data/monomicrobial_samples/Zymo-GridION-EVEN-BB-SN.fq.gz": "/output/assemblers/mock_microbial_community/",
    "/data/monomicrobial_samples/Zymo-GridION-LOG-BB-SN.fq.gz": "/output/assemblers/mock_microbial_community/"
}

def singularity_pull(tool, image):
    """
    Pull singularity images as docker images
    """
    print("-------------------------\n PULLING IMAGE FOR {}\n-------------------------------".format(tool))
    cmd = "module load apps/singularity && singularity pull {}.sif docker://{}".format(tool, image)
    process = subprocess.Popen(cmd, shell = True)
    process.wait()
    print("-----------------------")

def run_assembly(input_filepath, sample_name, threads, assembler, base_path, out_dir):
    """
    Run test datasets through assemblers with differing numbers of threads to assess scalability, performance and
    accuracy.
    """
    # scalability/performance
    print("--------------------------\n{} DE NOVO ASSEMBLY for {}\n--------------------------".format(assembler,
                                                                                                      input_filepath))
    os.chdir(out_dir)

    # create files for stderr and stdout
    stdout = out_dir + "stdout_{}_{}_{}_thread.txt".format(sample_name, assembler, threads)
    stderr = out_dir + "stderr_{}_{}_{}_thread.txt".format(sample_name, assembler, threads)

    # get command
    command = get_command(input_filepath, sample_name, assembler, threads, base_path, out_dir)

    if assembler == "raven":
        command += " > {}".format(stdout)

    # write command to file
    with open('command_file.txt', 'a') as command_file:
        command_file.writelines("SET OFF: " + command + "\n")
        command_file.close()

    if assembler == "raven":
        print(command)
        # wait for process to finish before printing returncode
        process = subprocess.Popen(command, shell=True, universal_newlines=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out, err = process.communicate()

    else:
        print(command)

        # wipe any existing file
        for file in stdout, stderr:
            outfile = open(file, 'w')
            outfile.close()

        with open(stdout, 'a') as out_file, open(stderr, 'a') as err_file:
            process = subprocess.Popen(command, shell = True, universal_newlines=True, bufsize=1,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # write stderr and stdout to file
            if stdout:
                for line in stdout:
                    sys.stdout.write(line)
                    out_file.write(line)
            if stderr:
                for line in stderr:
                    sys.stderr.write(line)
                    err_file.write(line)

            out_file.close()
            err_file.close()
        os.chdir(base_path)
    print("--------------------------")

def get_command(input_filepath, sample_name, assembler, threads, base_path, out_dir):
    """
    Returns command to run/append to file
    """
    assembly_prefix = "{}_{}_{}_thread".format(sample_name, assembler, threads)

    if 'Zymo' in input_filepath:
        genomeSize = "6.1944m"
    if 'enterococcus_faecium' in input_filepath:
        genomeSize = "2.85m"

    command_dictionary = {
        "flye": ("module load apps/singularity && /usr/bin/time -o {}time_{}.txt -v singularity exec "
                 "--bind `pwd`:`pwd` {}/{}.sif flye --nano-raw {} --out-dir {} --meta --threads "
                 "{}").format(out_dir, assembly_prefix, base_path, assembler, input_filepath, out_dir, threads),
        "canu": ("module load apps/singularity && /usr/bin/time -o {}time_{}.txt -v singularity exec "
                 "--bind `pwd`:`pwd` {}/{}.sif canu -p {} -d {} genomeSize={} maxThreads={} "
                 "maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 "
                 "-nanopore {}").format(out_dir, assembly_prefix, base_path, assembler, assembly_prefix, out_dir,
                                        genomeSize, threads, input_filepath),
        "raven": ("module load apps/singularity && /usr/bin/time -o {}time_{}.txt -v singularity exec "
                  "--bind `pwd`:`pwd` {}/{}.sif raven {} "
                  "-t {}").format(out_dir, assembly_prefix, base_path, assembler, input_filepath, threads),
        "canu_2_1": ("module load apps/singularity && /usr/bin/time -o {}time_{}.txt -v singularity exec "
                     "--bind `pwd`:`pwd` {}/{}.sif canu -p {} -d {} genomeSize={} maxThreads={} corThreads={} "
                     "maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 "
                     "-nanopore {}").format(out_dir, assembly_prefix, base_path, assembler, assembly_prefix, out_dir,
                                            genomeSize, threads, threads, input_filepath)
    }
    if assembler == "canu" and threads in (2,1):
        print(assembler, threads)
        return command_dictionary["canu_2_1"]
    else:
        return command_dictionary[assembler]

def main():
    base_path = os.getcwd().rsplit('/', 2)[0]
    # pull images
    for assembler in assembler_dictionary:
        singularity_pull(assembler, assembler_dictionary[assembler])

    # overwrite any existing file
    with open('command_file.txt', 'w') as command_file:
        command_file.close()

    # Runs assembly for each of the files for each number of threads for each assembler
    for assembler in assembler_dictionary:
        for fq in input_filepaths:
            for threads in [16, 8, 4, 2, 1]:
                # parse sample name
                sample_name = str(fq).rsplit("/", 1)[1].rsplit(".")[0]
                input_filepath = base_path + fq

                # create new directory for assembly if doesn't exist already
                out_dir = base_path + input_filepaths[fq] + "{}/{}/{}_{}/".format(sample_name, assembler, assembler,
                                                                                  threads)
                if os.path.isdir(out_dir):
                    print("Skipping: assembly already run")
                else:
                    os.makedirs(out_dir, exist_ok=True)
                    run_assembly(input_filepath, sample_name, threads, assembler, base_path, out_dir)

if __name__ == '__main__':
    main()

