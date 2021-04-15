import subprocess

file_path = "/mnt/flavia/rduffin/outbreak_pipeline/output/"
after_human_read_removal = ["200408_CMG_RUN3_Pool3/human_read_removal/200408_CMG_RUN3_Pool3_barcode01_unmapped.fastq",
                            "200416_CMG_RUN4_Pool6/human_read_removal/200416_CMG_RUN4_Pool6_barcode05_unmapped.fastq",
                            "200416_CMG_RUN4_Pool6/human_read_removal/200416_CMG_RUN4_Pool6_barcode08_unmapped.fastq",
                            "200504_CMG_RUN7_pool7/human_read_removal/200504_CMG_RUN7_pool7_barcode03_unmapped.fastq",
                            "200504_CMG_RUN7_pool7/human_read_removal/200504_CMG_RUN7_pool7_barcode06_unmapped.fastq",
                            "200611_pool11_CMG/human_read_removal/200611_pool11_CMG_barcode10_unmapped.fastq",
                            "200611_pool11_CMG/human_read_removal/200611_pool11_CMG_barcode11_unmapped.fastq",
                            "200617_pool13_cmg/human_read_removal/200617_pool13_cmg_barcode09_unmapped.fastq",
                            "200618_pool14_cmg/human_read_removal/200618_pool14_cmg_barcode01_unmapped.fastq",
                            "200618_pool14_cmg/human_read_removal/200618_pool14_cmg_barcode02_unmapped.fastq",
                            "200618_pool14_cmg/human_read_removal/200618_pool14_cmg_barcode04_unmapped.fastq",
                            "200618_pool14_cmg/human_read_removal/200618_pool14_cmg_barcode06_unmapped.fastq",
                            "200622_pool15_cmg/human_read_removal/200622_pool15_cmg_barcode08_unmapped.fastq",
                            "200622_pool15_cmg/human_read_removal/200622_pool15_cmg_barcode10_unmapped.fastq",
                            "200624_pool16_cmg/human_read_removal/200624_pool16_cmg_barcode10_unmapped.fastq",
                            "200729_POOL17/human_read_removal/200729_POOL17_barcode06_unmapped.fastq"]


def main():
    # overwrite old read count file
    read_count = open("read_count.txt", 'w')
    read_count.close()
    # get read count from all files in dictionary and print to output file
    with open("read_count.txt", 'a') as filetowrite:
        for file in after_human_read_removal:
            stdout, stderr = subprocess.Popen("grep '>' {}{} | wc -l ".format(file_path, file), shell=True,
                                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
            print(stdout.decode('ascii'))
            filetowrite.write("{}    {}".format(file, stdout.decode('ascii')))
    filetowrite.close()


if __name__ == '__main__':
    main()
