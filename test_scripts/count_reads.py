import subprocess
import os

file_dictionary = {"monoisolates": {"file_path": "/input/enterococcus_faecium/enterococcus/",
                                    "files": ["ef1_bc_75/210612_EF_R1_barcode01.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode02.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode03.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode04.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode05.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode06.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode07.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode08.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode09.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode10.fq",
                                              "ef1_bc_75/210612_EF_R1_barcode11.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode01.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode02.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode03.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode04.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode05.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode06.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode07.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode08.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode09.fq",
                                              "ef2_bc_75/210612_EF_R2_barcode10.fq"]
                                    },
                   "metagenomics": {"file_path": "/output/",
                                    "files": ["200408_CMG_RUN3_Pool3/human_read_removal/"
                                              "200408_CMG_RUN3_Pool3_barcode01_unmapped.fastq",
                                              "200416_CMG_RUN4_Pool6/human_read_removal/"
                                              "200416_CMG_RUN4_Pool6_barcode05_unmapped.fastq",
                                              "200416_CMG_RUN4_Pool6/human_read_removal/"
                                              "200416_CMG_RUN4_Pool6_barcode08_unmapped.fastq",
                                              "200504_CMG_RUN7_pool7/human_read_removal/"
                                              "200504_CMG_RUN7_pool7_barcode03_unmapped.fastq",
                                              "200504_CMG_RUN7_pool7/human_read_removal/"
                                              "200504_CMG_RUN7_pool7_barcode06_unmapped.fastq",
                                              "200611_pool11_CMG/human_read_removal/"
                                              "200611_pool11_CMG_barcode10_unmapped.fastq",
                                              "200611_pool11_CMG/human_read_removal/"
                                              "200611_pool11_CMG_barcode11_unmapped.fastq",
                                              "200617_pool13_cmg/human_read_removal/"
                                              "200617_pool13_cmg_barcode09_unmapped.fastq",
                                              "200618_pool14_cmg/human_read_removal/"
                                              "200618_pool14_cmg_barcode01_unmapped.fastq",
                                              "200618_pool14_cmg/human_read_removal/"
                                              "200618_pool14_cmg_barcode02_unmapped.fastq",
                                              "200618_pool14_cmg/human_read_removal/"
                                              "200618_pool14_cmg_barcode04_unmapped.fastq",
                                              "200618_pool14_cmg/human_read_removal/"
                                              "200618_pool14_cmg_barcode06_unmapped.fastq",
                                              "200622_pool15_cmg/human_read_removal/"
                                              "200622_pool15_cmg_barcode08_unmapped.fastq",
                                              "200622_pool15_cmg/human_read_removal/"
                                              "200622_pool15_cmg_barcode10_unmapped.fastq",
                                              "200624_pool16_cmg/human_read_removal/"
                                              "200624_pool16_cmg_barcode10_unmapped.fastq",
                                              "200729_POOL17/human_read_removal/"
                                              "200729_POOL17_barcode06_unmapped.fastq"]
                                    }
                   }


def count_reads(sample_type, fastq_list, file_path):
    """
    Calculate number of reads in each fastq file and output to file
    """
    print("COUNTING READS FOR SAMPLES: {}".format(sample_type))
    # get read count from all files in dictionary and print to output file
    with open("read_count_{}.txt".format(sample_type), 'w') as filetowrite:
        for file in fastq_list:
            command = "awk 'END {print NR/4}' " + file_path + file
            print(command)
            stdout, stderr = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT).communicate()
            filetowrite.write(stdout.decode('ascii'))

def calc_average_readno(sample_type):
    """
    Calculate average number of reads per sample for the sample type
    """
    print("CALCULATING AVERAGE READ COUNT FOR SAMPLES: {}".format(sample_type))

    # get read count from all files in dictionary and print to output file
    with open("average_reads_{}.txt".format(sample_type), 'w') as filetowrite:
        command = "cat read_count_{}.txt".format(sample_type) + "| jq -s add/length | awk '{x+=$0}END{print x/NR}'"
        print(command)
        stdout, stderr = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT).communicate()
        filetowrite.write(stdout.decode('ascii'))
    filetowrite.close()


def calculate_read_length(sample_type, fastq_list, file_path):
    """
    Calculate length of reads in each fastq file and output to file
    """
    print("CALCULATING READ LENGTHS FOR SAMPLES: {}".format(sample_type))
    # get max and min read length from all files in dictionary and print to output file
    with open("read_length_{}.txt".format(sample_type), 'w') as filetowrite:
        for file in fastq_list:
            command = "awk 'NR%4==2{print length($0)}' " + file_path + file
            print(command)
            stdout, stderr = subprocess.Popen(command, shell=True,
                                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
            filetowrite.write(stdout.decode('ascii'))
        filetowrite.close()


def calculate_length_summary(sample_type):
    """
    Output min, max, average and median read lengths and output to file.
    """
    print("GENERATING READ LENGTH SUMMARY FOR SAMPLES: {}".format(sample_type))

    with open("summary_read_length_{}.txt".format(sample_type), 'w') as filetowrite:
        command="cat read_length_" + sample_type + \
                ".txt|jq -s '{minimum:min,maximum:max,average:(add/length)," \
                "median:(sort|if length%2==1 then.[length/2|floor]else[.[length/2-1,length/2]]|add/2 end)," \
                "stdev:((add/length)as$a|map(pow(.-$a;2))|add/length|sqrt)}'"
        print(command)
        stdout, stderr = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT).communicate()

        filetowrite.write(stdout.decode('ascii'))
        filetowrite.close()

def main():
    for sample_type in file_dictionary:
        base_path = os.getcwd().rsplit('/', 2)[0]
        file_path = base_path + file_dictionary[sample_type]["file_path"]
        fastq_list = file_dictionary[sample_type]["files"]
        calculate_read_length(sample_type, fastq_list, file_path)
        calculate_length_summary(sample_type)
        count_reads(sample_type, fastq_list, file_path)
        calc_average_readno(sample_type)


if __name__ == '__main__':
    main()
