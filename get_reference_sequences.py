"""
Script to download required reference sequences from NCBI refseq for all scripts.
"""
import os
import gzip
import shutil
import config

def download_sequences(dictionary, out_dir, base_path):
    """
    Download all sequences in supplied dictionary.

    Checks if download already exists. If it doesn't, download the file, get the file name by splitting the dictionary
    value, and append bacterial name to start using the dictionary key.
    """
    print("--------------------------\nDOWNLOADING REFERENCE SEQUENCES\n--------------------------")
    for key in dictionary:
        file_name = str(dictionary[key]).rsplit("/", 1)[1]
        filepath= "{}{}_{}".format(out_dir, key, file_name)
        unzipped_file = str(filepath).rsplit(".gz", 1)[0]

        # download file
        if not os.path.exists(unzipped_file):
            print("Downloading: {} refseq file".format(key))
            command  = "wget -cO - {} > {}".format(dictionary[key], filepath)
            subprocess.run(command, shell=True)
            print("Download complete: {} refseq file".format(key))
            # unzip file
            with gzip.open(filepath, 'rb') as f_in:
                with open(unzipped_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # remove zipped file
            os.remove(filepath)
        else:
<<<<<<< HEAD
            print("REFERENCE FILE ALREADY DOWNLOADED FROM NCBI: {}".format(key))
=======
            print("REFERENCE FILE ALREADY DOWNLOADED FROM NCBI: {}".format(key))
        os.chdir(base_path)
    print("--------------------------")


def main():
    base_path = os.getcwd()
    out_dir = "{}/input/reference_sequences/".format(base_path)
    download_sequences(refseq_test_dict, out_dir, base_path)
    download_sequences(refseq_dict, out_dir, base_path)

if __name__ == '__main__':
    main()
>>>>>>> b782db751dd3e6341c28d65a3745f3ac3fed1774
