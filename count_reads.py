import subprocess

file_dictionary = {"monoisolates": {"file_path": "/mnt/flavia/rduffin/outbreak_pipeline/test_data/input/"
                                                 "enterococcus_faecium/enterococcus/",
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
                   "metagenomics": {"file_path": "/mnt/flavia/rduffin/outbreak_pipeline/output/",
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
    Calculate number of reads in each sample and output to file
    """
    # overwrite old read count file
    read_count = open("read_count_{}.txt".format(sample_type), 'w')
    read_count.close()
    # get read count from all files in dictionary and print to output file
    with open("read_count_{}.txt".format(sample_type), 'a') as filetowrite:
        for file in fastq_list:
            stdout, stderr = subprocess.Popen("grep '>' {}{} | wc -l ".format(file_path, file), shell=True,
                                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
            print(stdout.decode('ascii'))
            filetowrite.write("{}    {}".format(file, stdout.decode('ascii')))
    filetowrite.close()


def main():
    for key in file_dictionary:
        count_reads(sample_type=key, fastq_list=file_dictionary[key]["files"],
                    file_path=file_dictionary[key]["file_path"])


if __name__ == '__main__':
    main()
